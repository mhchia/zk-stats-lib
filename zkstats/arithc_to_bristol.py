from pathlib import Path
from dataclasses import dataclass
import json


# E.g. "gates": [
# { "op": "AMul", "lh_in": 6, "rh_in": 22, "out": 24 },
@dataclass(frozen=True)
class AGate:
  # just a serial id
  id: int
  # gate_type: AAdd, AMul, ...
  type: str
  # gate's lhs input wire id:
  lhs: int
  # gate's rhs input wire id:
  rhs: int
  # gate's output wire id:
  out: int


# E.g. "signals": {
#   "136": { "name": "TFMul.in[0][0]", "value": null },
# }
@dataclass(frozen=True)
class ASignal:
  signal_id: int
  name: str
  value: int


# "nodes": {
#   "762": { "is_const": false, "is_out": true, "signals": [496] },
# }
@dataclass(frozen=True)
class ANode:
  id: int
  # signal ids that this node is connected to
  signals: list[int]
  # is node a constant
  is_const: bool
  is_out: bool


# node.is_out=false
# [5]: "TFLog.e_until[0]"
# [0,3,51]
#   - 0: "0.in[0][0]"
#   - 3: "TFLog.in[0][0]"
#   - 51: "TFLog.x"
# [397]: "TFLog.x_over_b_minus_one_exp[0]"
#
# Plan: use the old way to determine
#   - inputs/outputs: old way. FIXME: need to figure out how `is_out` is determined
#   - constants: use signal.value != null


# signal["1"]: { "name": "0.tf_log_1_out[0]", "value": null },
# nodes["902"]: { "is_const": false, "is_out": true, "signals": [1, 621, 4, 2] },


def parse_arithc_json(
  arithc_path: Path,
  bristol_circuit_output_path: Path,
  circuit_info_output_path: Path,
):
  anodes, agates, main_inputs, main_outputs, const_names, const_values = _parse_arithc_json(arithc_path)
  _generate_bristol_and_circuit_info(bristol_circuit_output_path, circuit_info_output_path, agates, main_inputs, main_outputs, const_names, const_values)


def _parse_arithc_json(arithc_path: str):
  with open(arithc_path) as f:
    data = json.load(f)

  # signal_id -> ASignal
  asignals: dict[int, ASignal] = {}
  # node_id -> ANode
  anodes: dict[int, ANode] = {}
  # constant values: node_id -> const_value
  const_values: dict[int, int] = {}
  # constant names: node_id -> wire_name
  const_names: dict[int, str] = {}
  # highest level inputs: node_id -> wire_name
  main_inputs: dict[int, str] = {}
  # highest level outputs: node_id -> wire_name
  main_outputs: dict[int, str] = {}
  # gate_id -> AGate
  agates: dict[int, AGate] = {}

  for k, v in data['signals'].items():
    signal_id = int(k)
    signal_name = v['name']
    signal_value = int(v['value']) if v['value'] is not None else None
    asignal = ASignal(signal_id, signal_name, signal_value)
    asignals[signal_id] = asignal

  for k, v in data['nodes'].items():
    node_id = int(k)
    node_is_const = v['is_const']
    node_is_out = v['is_out']
    node_signal_ids = list(map(int, v['signals']))

    anode = ANode(id=node_id, signals=node_signal_ids, is_const=node_is_const, is_out=node_is_out)
    anodes[anode.id] = anode

  for gate_id, gate in enumerate(data['gates']):
    gate_type = gate['op']
    lh_input_node_id = gate['lh_in']
    rh_input_node_id = gate['rh_in']
    output_node_id = gate['out']
    agate = AGate(id=gate_id, type=gate_type, lhs=lh_input_node_id, rhs=rh_input_node_id, out=output_node_id)
    agates[gate_id] = agate

  # find
  # 1. highest level inputs
  #   - signals name starts with "0."
  # 2. highest level outputs
  # 3. constant values for each node

  for node_id, anode in anodes.items():
    # FIXME: this seems not a correct way to determine input/output
    # Now we assume if the node name starts with "0.", it's a highest level signal.
    # If there are other information in the future, we should use that instead.
    signal_ids = anode.signals
    node_signals = [asignals[sid] for sid in signal_ids]

    # Iterate through all signals of this node to find the highest level inputs and outputs
    for _signal in node_signals:
      signal_name = _signal.name
      if signal_name.startswith("0."):
        # Assuming inputs are always before outputs
        # Only set it if it's not already set for inputs, to avoid names being overwritten
        if anode.id not in main_inputs:
          main_inputs[anode.id] = signal_name[2:]
        # It's fine for output to be overwritten
        main_outputs[anode.id] = signal_name[2:]

    # Iterate through all signals of this node to find the constant values
    # FIXME: are these sanity checks correct?
    # Sanity check 1: if a node is a constant, all its signals must have value != null
    if anode.is_const:
      is_one_signal_has_value = any(_signal.value is not None for _signal in node_signals)
      if not is_one_signal_has_value:
        raise Exception(f"Constant node has no signal with value not None: {anode=}")
    # Sanity check 2: there can be only one signal with value != null
    is_multiple_signals_have_value = sum(1 for _signal in node_signals if _signal.value is not None) > 1
    if is_multiple_signals_have_value:
      raise Exception(f"Node has multiple signals with value not None: {anode=}")
    # If any of the node's signals has value, set it as the node's constant value
    for _signal in node_signals:
      if _signal.value is not None:
        const_values[anode.id] = _signal.value
        const_names[anode.id] = signal_name

  # Clean up `main_inputs` and `main_outputs` to make sure they only contain their corresponding wires
  # gid = 0...ngate-1
  # A node can be simultaneously an input and an output
  # If it's an input to a gate, remove it from `main_outputs`
  for gid in agates:
    gate = agates[gid]
    lhs = gate.lhs
    rhs = gate.rhs
    out = gate.out
    # Remove lhs input from `main_outputs`
    main_outputs.pop(lhs, None)
    # Remove rhs input from `main_outputs`
    main_outputs.pop(rhs, None)
    # Remove output from `main_inputs`
    main_inputs.pop(out, None)

  # Now we get the correct `main_inputs` and `main_outputs`
  print("!@# after pop: main_inputs= ", main_inputs)
  print("!@# after pop: main_outputs=", main_outputs)
  const_name_to_values = {const_names[node_id]: const_value for node_id, const_value in const_values.items()}
  print("!@# after pop: const name to values=  ", const_name_to_values)
  return anodes, agates, main_inputs, main_outputs, const_names, const_values


def _generate_bristol_and_circuit_info(
  bristol_circuit_path: str,
  circuit_info_filepath: str,
  agates: dict[int, AGate],
  main_inputs: dict[int, str],
  main_outputs: dict[int, str],
  const_names: dict[int, str],
  const_values: dict[int, int],
) -> tuple[str, str]:
  '''
  To assign nodes (wires) for bristol circuit, we need to make each input wire has a lower index than its output
  To achieve this, we can assign index for each wire in their topological order:
    1. Based on the parsed arithc, build a tree structure for the this circuit
      - Each node (wire) is a tree node, output is the parent of its inputs
      - There can be multiple roots (outputs)
    2. Traverse the tree in topological order and assign index for each wire
  '''
  class TNode:
    def __init__(self, rid: int, type: str | None, lnode: 'TNode', rnode: 'TNode'):
      # id in the topological order (wire id in bristol circuit)
      self.iid = 0
      # original node id from arithc
      self.rid = rid
      self.lnode = lnode
      self.rnode = rnode
      # gate type: AAdd, AMul. None if it's an input or output without a gate
      self.type = type
      self.is_root = True
      self.is_leaf = True

  class TTree:
    def __init__(self):
      self.roots: list[TNode] = []
      # node_id -> TNode. All nodes wires
      self.tnodes: dict[int, TNode] = {}
      self.leaves: dict[int, TNode] = {}
      self.gcount = 0
      self.wcount = 0

    def insert_node(self, gate: AGate):
      gtype = gate.type
      gout = gate.out
      # left wire
      glhs = gate.lhs
      # right wire
      grhs = gate.rhs
      # TNode(rid, type, lnode, rnode)
      tnode = TNode(gout, gtype, None, None)
      self.tnodes[gout] = tnode
      # if left wire is not in the tree, create a new node for it
      if (self.tnodes.get(glhs) is None):
        # wire
        lnode = TNode(glhs, None, None, None)
        self.tnodes[glhs] = lnode
      # if right wire is not in the tree, create a new node for it
      if (self.tnodes.get(grhs) is None):
        rnode = TNode(grhs, None, None, None)
        self.tnodes[grhs] = rnode

    def build_tree(self, gates: dict[int, AGate]):
      # `gates` is ordered by id.
      # Insert all wires into self.tnodes
      for gid in gates:
        gate = gates[gid]
        self.insert_node(gate)

      # for each gate, set input wires as children of the output wire
      # also, count gates `self.gcount` and wires `self.wcount`
      for gid in gates:
        gate = gates[gid]
        gout = gate.out
        # output node
        tnode = self.tnodes[gout]
        glhs = gate.lhs
        grhs = gate.rhs
        lnode = self.tnodes[glhs]
        rnode = self.tnodes[grhs]
        lnode.is_root = False
        rnode.is_root = False
        tnode.is_leaf = False
        tnode.lnode = lnode
        tnode.rnode = rnode
        self.gcount += 1

      # If a node is a root, add it to `self.roots`
      # If a node is a leaf, add it to `self.leaves`
      for tnid in self.tnodes:
        tnode = self.tnodes[tnid]
        assert tnode.rid == tnid
        if tnode.is_root:
          self.roots.append(tnode)
        elif tnode.is_leaf:
          self.leaves[tnid] = tnode

      self.topological_sort_nodes()
      self.wcount = len(self.sorted_wires)

    def topological_sort_nodes(self):
      # Count how many wires a node is pointed to
      pointed_to = {}
      for nid in self.tnodes:
        pointed_to[nid] = 0
      for nid in self.tnodes:
        tnode = self.tnodes[nid]
        if tnode.lnode:
          pointed_to[tnode.lnode.rid] += 1
        if tnode.rnode:
          pointed_to[tnode.rnode.rid] += 1
      # print("!@# before: pointed_to=", pointed_to)
      # Start from the known roots (outputs)
      reversed_roots = self.roots[::-1]
      for i in reversed_roots:
        pointed_to.pop(i.rid)
      # print("!@# after:  pointed_to=", pointed_to)
      queue = list(reversed_roots)
      # print("!@# queue before while=", queue)
      wires: list[TNode] = []
      wire_index = 0
      while queue:
        # Visit next node
        current_node = queue.pop(0)
        # print("")
        # print("")
        # print("!@# loop i=", wire_index, ": visit current_node.rid=", current_node.rid)
        # print("!@# loop i=", wire_index, ": queue before decrement=", [i.rid for i in queue])
        # print("!@# loop i=", wire_index, ": pointed_to before decrement=", pointed_to)
        current_node.iid = wire_index
        wires.append(current_node)
        # Decrement the count of the nodes it points to
        if current_node.lnode is not None:
          if pointed_to[current_node.lnode.rid] <= 0:
            raise Exception("pointed_to <= 0")
          pointed_to[current_node.lnode.rid] -= 1
          # If the count is 0, add it to the queue and delete it from `pointed_to`
          if pointed_to[current_node.lnode.rid] == 0:
            queue.append(current_node.lnode)
            pointed_to.pop(current_node.lnode.rid)
          elif pointed_to[current_node.lnode.rid] < 0:
            raise Exception("pointed_to < 0")
        if current_node.rnode is not None:
          if pointed_to[current_node.rnode.rid] <= 0:
            raise Exception("pointed_to <= 0")
          pointed_to[current_node.rnode.rid] -= 1
          # If the count is 0, add it to the queue and delete it from `pointed_to`
          if pointed_to[current_node.rnode.rid] == 0:
            queue.append(current_node.rnode)
            pointed_to.pop(current_node.rnode.rid)
          elif pointed_to[current_node.rnode.rid] < 0:
            raise Exception("pointed_to < 0")
        # print("!@# loop i=", wire_index, ": queue after decrement=", [i.rid for i in queue])
        # print("!@# loop i=", wire_index, ": pointed_to after decrement=", pointed_to)
        wire_index += 1
      # Now wires contains all nodes in topological order, i.e. every parent (output) has lower index
      # than their children (inputs). Since we want inputs having lower values than outputs, reverse the list
      self.sorted_wires = wires[::-1]
      for index, node in enumerate(self.sorted_wires):
        node.iid = index
      # self.map_node_rid_to_wire_index = {i.rid: index for index, i in enumerate(self.sorted_wires)}

    def generate_bristol_circuit(self):
      num_gates = self.gcount
      num_wires = self.wcount
      num_inputs = len(self.leaves)
      num_outputs = len(self.roots)
      header = f"""{num_gates} {num_wires}
{num_inputs} {' '.join(["1"]* num_inputs)}
{num_outputs} {' '.join(["1"]* num_outputs)}

"""
      # Each gate takes a line: "2 1 1 2 3 AAdd"
      gates = "\n".join([
        f"2 1 {node.lnode.iid} {node.rnode.iid} {node.iid} {node.type}"
        for node in self.sorted_wires
        if node.type is not None
      ])
      return header + gates

    def generate_circuit_info(self):
      #
      # Write circuit_info.json:
      # {
      #   "input_name_to_wire_index": ,
      #   "constants":,
      #   "output_name_to_wire_index",
      # }

      rid_to_iid = {node.rid: node.iid for node in tt.sorted_wires}
      # Map input name to wire index in MP-SPDZ circuit (including constant wires)
      input_name_to_wire_index = {
        input_name: rid_to_iid[node_rid]
        # for node_rid in self.leaves if node_rid not in const_values
        for node_rid, input_name in main_inputs.items()
      }

      # FIXME: outputs without a gate are skipped (i.e. direct assigned from input or a constant, etc)

      # Prepare constants: const_values is what we want
      # Just sanity check for all constant must be in leaves so we don't miss passing any of them to MP-SPDZ circuit
      const_name_to_value_wire_id = {
        const_names[node_rid]: {
          'value': const_value,
          'wire_index': rid_to_iid[node_rid],
        }
        for node_rid, const_value in const_values.items()
        if node_rid in rid_to_iid  # Skip constant wires that are not used in any gates. E.g. constant outputs
      }

      # Prepare outputs
      # Map output name to wire index in MP-SPDZ circuit
      output_name_to_wire_index = {
        output_name: rid_to_iid[node_rid]
        for node_rid, output_name in main_outputs.items()
        if node_rid in rid_to_iid  # Skip output wires that are not used in any gates. E.g. constant outputs
      }
      print("!@# output_name_to_wire_index=", output_name_to_wire_index)
      return {
        "input_name_to_wire_index": input_name_to_wire_index,
        "constants": const_name_to_value_wire_id,
        "output_name_to_wire_index": output_name_to_wire_index,
      }

  tt = TTree()
  tt.build_tree(agates)

  with open(bristol_circuit_path, "w") as f:
    f.write(tt.generate_bristol_circuit())

  with open(circuit_info_filepath, "w") as f:
    json.dump(tt.generate_circuit_info(), f, indent=4)


def main():
  project_root = Path(__file__).resolve().parent
  arithc_path = project_root / "output-path" / "circuit.json"

  with open(arithc_path) as f:
    data = json.load(f)

  # skip `node_count` for now. It should subtract the number of removed nodes?
  try:
    node_count = data['node_count']
    print("!@# node_count=", node_count)
  except KeyError:
    pass

  try:
    _vars = data['vars']
    print("!@# len(vars)=", len(_vars))
  except KeyError:
    pass
  try:
    signals = data['signals']
    print("!@# len(signals)=", len(signals))
  except KeyError:
    pass

  try:
    nodes = data['nodes']
    print("!@# len(nodes)=", len(nodes))
  except KeyError:
    pass

  try:
    gates = data['gates']
    print("!@# len(gates)=", len(gates))
  except KeyError:
    pass

  bristol_circuit_output_path = arithc_path.parent / f"{arithc_path.stem}.txt"
  circuit_info_filepath = arithc_path.parent / f"{arithc_path.stem}.circuit_info.json"

  parse_arithc_json(arithc_path, bristol_circuit_output_path, circuit_info_filepath)

if __name__ == '__main__':
  main()
