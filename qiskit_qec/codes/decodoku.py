from random import choice
from numpy.random import poisson

from rustworkx.visualization import mpl_draw

from qiskit_qec.utils.visualizations import QiskitGameEngine
from qiskit_qec.decoders import DecodingGraph


class Decodoku:
    def __init__(self, p=0.1, k=2, d=10, process=None, errors=None, nonabelian=False):

        self.p = p
        self.k = k
        self.d = d
        self.process = process
        self.decoder = self.process != None
        if errors:
            self.errors = errors
        else:
            self.errors = []
        self.nonabelian = nonabelian

        self.L = d + 1

        self.node_color = []
        self.boundary_corrections = []

        self.generate_syndrome()
        self.generate_graph()

    def generate_syndrome(self):

        syndrome = {}
        for x in range(self.L):
            for y in range(self.L):
                syndrome[x, y] = 0

        if self.errors:
            error_num = len(self.errors)
        else:
            error_num = poisson(self.p * 2 * self.L**2)
            for _ in range(error_num):

                x0 = choice(range(self.L))
                y0 = choice(range(self.L))

                if choice([True, False]):
                    dx = []
                    if x0 > 0:
                        dx.append(-1)
                    if x0 < (self.L - 1):
                        dx.append(+1)
                    x1 = x0 + choice(dx)
                    y1 = y0
                else:
                    dy = []
                    if y0 > 0:
                        dy.append(-1)
                    if y0 < (self.L - 1):
                        dy.append(+1)
                    x1 = x0
                    y1 = y0 + choice(dy)
                self.errors.append((x0, y0, x1, y1))

        for x0, y0, x1, y1 in self.errors:
            if x0 not in [0, self.L - 1] or x1 not in [0, self.L - 1]:
                e = choice(range(1, self.k))
                syndrome[x0, y0] += e
                syndrome[x1, y1] += self.k - e

        # generate colours for clusters
        self.error_colors = [
            "#00ff00",
            "#fc04c5",
            "#00ffff",
            "#007f7f",
            "#bc020f",
            "#767202",
            "#8156f8",
            "#51ee7c",
            "#ee667d",
            "#55017c",
            "#007fff",
            "#dcfe44",
            "#ff7fff",
            "#7f7f7f",
            "#7fff00",
            "#5fb6df",
            "#009f08",
            "#7f00ff",
            "#fe4b0b",
            "#163747",
            "#0330aa",
            "#7fffff",
            "#a6b231",
            "#a4369c",
            "#eeae74",
            "#af4c3c",
            "#01c6c1",
            "#692801",
            "#fa2261",
            "#02da4d",
            "#52b846",
            "#db45e1",
            "#1c3ef8",
            "#3c72bd",
            "#055902",
            "#a58eca",
            "#99c689",
            "#41f4ca",
            "#b2006a",
            "#387144",
            "#6d365e",
            "#f9ef98",
            "#00007f",
            "#410dc2",
            "#3cb490",
            "#b1f9a7",
            "#8dec4c",
            "#d87c1c",
            "#03fa97",
            "#4dce00",
            "#ff0000",
            "#ffff00",
            "#c113f0",
            "#c7d808",
            "#8e02b1",
            "#3c0538",
            "#633eb1",
            "#1ec5fd",
            "#e790b8",
            "#38fe31",
            "#a8b6fd",
            "#a47e47",
            "#2a4d7f",
            "#009844",
            "#527bfd",
            "#840336",
            "#fdb7ee",
            "#c43304",
            "#c8d45f",
            "#7ddec0",
            "#f63b9f",
            "#c2f8f9",
            "#232c0b",
            "#b8968b",
            "#77a405",
            "#5f29f5",
            "#4f452a",
            "#0890b2",
            "#b764ac",
            "#3e8c0a",
            "#f4d431",
            "#0259cd",
            "#ce3c74",
            "#09ce12",
            "#6a9aac",
            "#962fd7",
            "#cc05a1",
            "#407f7a",
            "#dda737",
            "#89d204",
            "#b06efa",
            "#6b8836",
            "#04c380",
            "#21a0dc",
            "#7e66bf",
            "#a05505",
            "#84a763",
            "#fd1ffe",
            "#fb5149",
            "#c02140",
            "#f98148",
            "#7f0000",
            "#010f34",
            "#d197f2",
            "#c1a200",
            "#430501",
            "#e7033e",
            "#8194fe",
            "#53ddfc",
            "#a1597a",
            "#06f8c6",
            "#f7c3aa",
            "#0e08ac",
            "#7f6037",
            "#087427",
            "#252373",
            "#81fd8f",
            "#93b7c7",
            "#005d54",
            "#3e5c01",
            "#68da2e",
            "#fc67bb",
            "#ff007f",
            "#f5e5d2",
            "#635a76",
            "#4402f3",
            "#fe2422",
            "#37d0b4",
            "#353ccc",
            "#2bac25",
            "#1bfe6a",
            "#8e2c2c",
            "#9c2762",
            "#f9d569",
            "#ff7f00",
            "#081ad7",
            "#2ca765",
            "#83d0f9",
            "#db74dc",
            "#d16b52",
            "#392d9e",
            "#34d331",
            "#7a1f92",
            "#6ac1a1",
            "#b1fc0f",
            "#c8fc76",
            "#c8d998",
            "#aa7b01",
            "#d225c7",
            "#6f85d6",
            "#34cd6e",
            "#a8e7d7",
            "#bdadb4",
            "#7311d5",
            "#bfad63",
            "#7ccd63",
            "#1062a2",
            "#2965fa",
            "#d3521e",
            "#a652d8",
            "#63fd56",
            "#504dfd",
            "#fb8987",
            "#2ff4fc",
            "#d8ccfe",
            "#023f20",
            "#53fb0d",
            "#3897ad",
            "#651f3d",
            "#04fa34",
            "#79b82b",
            "#36f7a0",
            "#b63cf9",
            "#7d0a69",
            "#147ad5",
            "#a7d731",
            "#81519a",
            "#cf47a2",
            "#459648",
            "#63a481",
            "#fc59f0",
            "#62fbad",
            "#e33240",
            "#d5906a",
            "#2aea03",
            "#443251",
            "#21055f",
            "#cff6c9",
            "#dcbb11",
            "#5d6f9e",
            "#579ffd",
            "#565dd8",
            "#8afebd",
            "#aa7b73",
            "#f8d802",
            "#047c00",
            "#47e151",
            "#d4166d",
            "#14b3a5",
            "#949221",
        ]
        if len(self.error_colors) < error_num:
            self.error_colors = self.error_colors * (1 + int(error_num / len(self.error_colors)))

        # compute boundary parities and scrub their syndromes
        parity = [0, 0]
        for e, x in enumerate([0, self.L - 1]):
            for y in range(self.L):
                parity[e] += syndrome[x, y]
                syndrome[x, y] = 0
            parity[e] = parity[e] % self.k
        self.boundary_errors = parity

        self.syndrome = syndrome
        self.original_syndrome = syndrome.copy()
        self.boundary_errors = parity

    def generate_graph(self):

        dg = DecodingGraph(None)

        d = self.L - 1

        pos = []
        for x in range(1, self.L - 1):
            for y in range(self.L):
                e = x - 1
                t = self.L - 1 - y
                dg.graph.add_node({"y": t, "x": e, "is_boundary": False})
                pos.append((e, -t))
        for e in [0, 1]:
            dg.graph.add_node(
                {"y": 0, "x": (d - 1) * (e == 1) - 1 * (e == 0), "element": e, "is_boundary": True}
            )
            pos.append((d * (e == 1) - 2 * (e == 0), -(self.L - 1) / 2))

        nodes = dg.graph.nodes()
        # connect edges to boundary nodes
        for y in range(self.L):
            t = y
            n0 = nodes.index({"y": 0, "x": -1, "element": 0, "is_boundary": True})
            n1 = nodes.index({"y": t, "x": 0, "is_boundary": False})
            dg.graph.add_edge(n0, n1, None)
            n0 = nodes.index({"y": 0, "x": d - 1, "element": 1, "is_boundary": True})
            n1 = nodes.index({"y": t, "x": d - 2, "is_boundary": False})
            dg.graph.add_edge(n0, n1, None)
        # connect bulk nodes with space-like edges
        for y in range(self.L):
            for x in range(1, self.L - 2):
                t = y
                e = x - 1
                n0 = nodes.index({"y": t, "x": e, "is_boundary": False})
                n1 = nodes.index({"y": t, "x": e + 1, "is_boundary": False})
                dg.graph.add_edge(n0, n1, None)
        # connect bulk nodes with time-like edges
        for y in range(self.L - 1):
            for x in range(1, self.L - 1):
                t = y
                e = x - 1
                n0 = nodes.index({"y": t, "x": e, "is_boundary": False})
                n1 = nodes.index({"y": t + 1, "x": e, "is_boundary": False})
                dg.graph.add_edge(n0, n1, None)

        self.decoding_graph = dg
        self.graph_pos = pos
        self.update_graph()

    def update_graph(self, original=False):

        for node in self.decoding_graph.graph.nodes():
            node["highlighted"] = False
            if self.k != 2:
                node["value"] = 0

        if original:
            syndrome = self.original_syndrome
        else:
            syndrome = self.syndrome

        highlighted_color = []
        for node in self.decoding_graph.graph.nodes():
            if node["is_boundary"]:
                highlighted_color.append("orange")
            else:
                x = node["x"] + 1
                y = node["y"]
                if (syndrome[x, y] % self.k) > 0:
                    highlighted_color.append("red")
                    node["value"] = syndrome[x, y] % self.k
                    node["highlighted"] = True
                else:
                    if self.decoder:
                        highlighted_color.append("#bbbbbb")
                    else:
                        highlighted_color.append("cornflowerblue")
        self.node_color = highlighted_color

    def start(self, engine):

        d = self.k
        syndrome = self.syndrome

        # set edges to orange
        for x in [0, engine.L - 1]:
            for y in range(engine.L):
                engine.screen.pixel[x, y].set_text("")
                engine.screen.pixel[x, y].set_color("orange")

        # set bulk to blue
        for x in range(1, engine.L - 1):
            for y in range(engine.L):
                engine.screen.pixel[x, y].set_text("")
                engine.screen.pixel[x, y].set_color("blue")

        # display changed syndromes
        for x in range(1, engine.L - 1):
            for y in range(engine.L):
                if syndrome[x, y] % d != 0:
                    if x not in [0, engine.L - 1]:
                        if d != 2:
                            if not self.nonabelian:
                                engine.screen.pixel[x, y].set_text(str(syndrome[x, y] % d))
                            else:
                                if syndrome[x, y] % d == d / 2:
                                    engine.screen.pixel[x, y].set_text("Λ")
                                else:
                                    engine.screen.pixel[x, y].set_text("Φ")

                        engine.screen.pixel[x, y].set_color("red")

        engine.screen.pixel["text"].set_text("Choose node with an error")

    # this is the function that does everything
    def next_frame(self, engine):

        d = self.k
        syndrome = self.syndrome

        if len(engine.pressed_pixels) == 1:
            (x, y) = engine.pressed_pixels[0]
            if x not in [0, engine.L - 1] and syndrome[x, y] > 0:
                engine.screen.pixel["text"].set_text("Choose node to move it to")
            else:
                engine.pressed_pixels = []
        elif len(engine.pressed_pixels) == 2:
            # if we have a node and somewhere to move it to, then move it
            [(x0, y0), (x1, y1)] = engine.pressed_pixels
            if (x0, y0) != (x1, y1):
                syndrome[x1, y1] += syndrome[x0, y0]
                syndrome[x0, y0] = 0
                for x, y in [(x0, y0), (x1, y1)]:
                    if syndrome[x, y] % d == 0:
                        if x not in [0, engine.L - 1]:
                            engine.screen.pixel[x, y].set_text("")
                            engine.screen.pixel[x, y].set_color("blue")
                    else:
                        if x not in [0, engine.L - 1]:
                            if d != 2:
                                if not self.nonabelian:
                                    engine.screen.pixel[x, y].set_text(str(syndrome[x, y] % d))
                                else:
                                    if syndrome[x, y] % d == d / 2:
                                        engine.screen.pixel[x, y].set_text("Λ")
                                    else:
                                        engine.screen.pixel[x, y].set_text("Φ")
                            engine.screen.pixel[x, y].set_color("red")

                engine.pressed_pixels = []
                engine.screen.pixel["text"].set_text("Choose syndrome element")
            else:
                engine.pressed_pixels = [(x0, y0)]

        else:
            engine.pressed_pixels = []
            engine.screen.pixel["text"].set_text("Choose syndrome element")

        if engine.controller["next"].value:
            engine.pressed_pixels = []
            engine.screen.pixel["text"].set_text("Choose syndrome element")

        # see how many non-trivial syndromes are left
        num_elems = 0
        for x in range(1, engine.L - 1):
            for y in range(engine.L):
                num_elems += syndrome[x, y] % d

        if num_elems == 0:
            parity = [0, 0]
            for e, x in enumerate([0, engine.L - 1]):
                for y in range(engine.L):
                    parity[e] += syndrome[x, y]
                parity[e] = parity[e] % d
                self.boundary_corrections = parity
            for e, x in enumerate([0, engine.L - 1]):
                engine.screen.pixel[x, 0].set_text(str(self.boundary_errors[e]))
                engine.screen.pixel[x, self.L - 1].set_text(str(self.boundary_corrections[e]))

            if (self.boundary_corrections[0] + self.boundary_errors[0]) % d == 0:
                engine.screen.pixel["text"].set_text("Correction successful! :D")
                for y in range(engine.L):
                    for x in [0, engine.L - 1]:
                        engine.screen.pixel[x, y].set_color("green")
            else:
                engine.screen.pixel["text"].set_text("Correction unsuccessful! :(")
                for y in range(engine.L):
                    for x in [0, engine.L - 1]:
                        engine.screen.pixel[x, y].set_color("red")

            if engine.controller["next"].value:
                self.restart()
                engine.start(engine)

        self.update_graph()

    def restart(self):
        self.errors = []
        self.generate_syndrome()
        self.update_graph()

    def draw_graph(self, clusters=True):
        if self.decoder and clusters:
            parity, clusters = self.process(self)
        else:
            parity, clusters = None, None

        node_color = self.node_color.copy()
        if clusters:
            for n, c in clusters.items():
                node_color[n] = self.error_colors[c]

        def get_label(node):
            if node["is_boundary"] and parity:
                return str(parity[node["element"]])
            elif node["highlighted"] and "value" in node and self.k != 2:
                return str(node["value"])
            else:
                return ""

        return mpl_draw(
            self.decoding_graph.graph,
            pos=self.graph_pos,
            node_color=node_color,
            labels=get_label,
            with_labels=True,
        )

    def run(self):
        return QiskitGameEngine(self.start, self.next_frame, L=self.L, active_screen=True)
