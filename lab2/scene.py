from manim import *
from manim.opengl import *  # type: ignore


class Bresenham(Scene):
    def prepare_grid(self):
        print('Preparing to set up grid')
        self.grid = []
        self.grid_mobj = VGroup()
        self.grid_coords = dict()
        square_size = 0.1

        print('Measuring screen')
        screen = FullScreenRectangle()

        height = screen.height
        width = screen.width
        x = 0
        y = 0

        maxx = 0
        maxy = 0
        square = Square(square_size, stroke_color=GRAY)
        square.stroke_width = 0.005
        square.align_to(screen, UP).align_to(screen, LEFT)
        for xr in range(142):
            print('xr', xr)
            for yr in range(80):
                s = square.copy()
                self.grid.append(s)
                self.grid_coords[(x, y)] = s
                self.grid_mobj.add(s)

                y += 1
                maxx = max(maxx, x)
                maxy = max(maxy, y)
                square.shift(square_size * DOWN)

            square.shift((height) * UP + square_size * RIGHT)
            y = 0
            x += 1
        
        print(maxx,maxy)

        self.grid_mobj.center()
        #self.play(LaggedStartMap(Write, self.grid, run_time=1))

        print('Making optimized grid')
        optimized_grid = VGroup()
        grid_xstep=square_size
        grid_ystep=square_size


        v = screen.get_vertices()[1]
        grid_xstep = abs(grid_xstep)
        count = int(self.grid_mobj.width / grid_xstep)
        grid = VGroup(
            *(
                Line(
                    v + i * grid_xstep * RIGHT,
                    v + i * grid_xstep * RIGHT + height * DOWN,
                    stroke_width=1,
                    color=GRAY,
                    stroke_opacity=0.5
                )
                for i in range(1, count)
            )
        )
        optimized_grid.add(grid)

        grid_ystep = abs(grid_ystep)
        count = int(self.grid_mobj.height / grid_ystep)
        grid = VGroup(
            *(
                Line(
                    v + i * grid_ystep * DOWN,
                    v + i * grid_ystep * DOWN + width * RIGHT,
                    color=GRAY,
                    stroke_width=1,
                    stroke_opacity=0.5
                )
                for i in range(1, count)
            )
        )
        optimized_grid.add(grid)

        print('Animating...')
        self.play(Create(optimized_grid))
        

    def turn_on(self, sq: Square, fast=False):        
            return AnimationGroup(
                Flash(sq, run_time=0.4 if not fast else 0.1),
                Succession(
                    # sq.animate(rate_func=rate_functions.rush_into, run_time=0.2).set_fill(color.YELLOW, opacity=1)),
                    sq.animate(rate_func=rate_functions.rush_from, run_time=0.4 if not fast else 0.1).set_fill(color.YELLOW_E, opacity=0.5)
                )
            )

        

    def draw_line(self, x0, y0, x1, y1, fast=False):
        speed = 1 if fast else 0.5
        sq1: Square = self.grid_coords[(x0, y0)]
        sq2: Square = self.grid_coords[(x1, y1)]
        ideal_line = Line(sq1.get_center(), sq2.get_center())
        self.play(Create(ideal_line), run_time=speed)

        said = []
        displayed = [ideal_line]
        def say(text):
            if fast: return
            tex = MathTex(text)
            if not said:
                tex.move_to(FullScreenRectangle().get_boundary_point(UP+LEFT)).shift_onto_screen()
            else:
                tex.next_to(said[-1], DOWN).align_to(said[-1], LEFT)
            said.append(tex)
            self.play(FadeIn(tex, shift=UP, run_time=0.5))
        
        def reset_say():
            if fast: return
            self.play(LaggedStartMap(lambda x: FadeOut(x, shift=LEFT), said, run_time=0.3))   
            said.clear() 


        try:
            say(f"A = ({x0}, {y0})")
            self.play(self.turn_on(sq1, fast))
            say(f"B = ({x1}, {y1})")

            dx = abs(x0 - x1)
            dy = abs(y0 - y1)
            say(f"dx = {dx}, dy = {dy}")

            if dx == 0:
                say('\\rightarrow vertical')
                a,b = sorted([y0, y1])
                for i in range(a, b+1):
                    self.play(self.turn_on(self.grid_coords[(x0, i)], fast))
                return
            
            say("S = \\frac{dy}{dx} =" + str(round(dy / dx, 3)))

            def linehigh(x0,y0,x1,y1):
                sq1: Square = self.grid_coords[(x0, y0)]
                sq2: Square = self.grid_coords[(x1, y1)]
                dx = x1-x0
                dy = y1-y0
                xi = 1
                xis = '+1'
                if dx < 0:
                    xi = -1
                    xis = '-1'
                    dx = -dx

                error_label = Tex("D = ").shift(3*DOWN)
                error_value = Integer(2*dx - dy).next_to(error_label)
                error = VGroup(error_label, error_value)
                displayed.append(error)
                self.play(Write(error), run_time=speed)

                selector = sq1.copy().set_fill(RED, 0.6)
                displayed.append(selector)
                self.play(Create(selector), run_time=speed)

                sel_x = Integer(x0)
                lx = Tex("x = ")
                sel_y = Integer(y0)
                ly = Tex("y = ")
                sel_coord = VGroup(lx, sel_x, ly, sel_y).arrange().shift(3*UP)
                displayed.append(sel_coord)
                self.play(Write(sel_coord), run_time=speed)

                sels = selector.width

                while sel_y.get_value() != y1:
                    reset_say()
                    self.play(self.turn_on(self.grid_coords[(sel_x.get_value(),sel_y.get_value())], fast))
                    D = error_value.get_value()
                    if D>0:
                        say(f"D > 0 \\rightarrow x = x {xis}")
                        self.play(AnimationGroup(
                            error_value.animate.set_value(D + (2*(dx - dy))),
                            sel_x.animate.set_value(sel_x.get_value() + xi),
                            sel_y.animate.set_value(sel_y.get_value() + 1),
                            selector.animate.shift((sels * RIGHT * xi) + (sels * DOWN))
                        ), run_time=speed)
                    else:
                        say("D \\leq 0 \\rightarrow x=x")
                        self.play(AnimationGroup(
                            error_value.animate.set_value(D + (2*dx)),
                            sel_y.animate.set_value(sel_y.get_value() + 1),
                            selector.animate.shift(sels * DOWN)
                        ), run_time=speed)
                self.play(self.turn_on(self.grid_coords[(sel_x.get_value(),sel_y.get_value())], fast))


            def linelow(x0,y0,x1,y1):
                sq1: Square = self.grid_coords[(x0, y0)]
                sq2: Square = self.grid_coords[(x1, y1)]
                dx = x1-x0
                dy = y1-y0
                yi = 1
                yis = '+1'
                if dy < 0:
                    yi = -1
                    yis = '-1'
                    dy = -dy

                error_label = Tex("D = ").shift(3*DOWN)
                error_value = Integer((2*dy) - dx).next_to(error_label)
                error = VGroup(error_label, error_value)
                displayed.append(error)
                self.play(Write(error), run_time=speed)

                selector = sq1.copy().set_fill(RED, 0.6)
                displayed.append(selector)
                self.play(Create(selector), run_time=speed)

                sel_x = Integer(x0)
                lx = Tex("x = ")
                sel_y = Integer(y0)
                ly = Tex("y = ")
                sel_coord = VGroup(lx, sel_x, ly, sel_y).arrange().shift(3*UP)
                displayed.append(sel_coord)
                self.play(Write(sel_coord), run_time=speed)

                sels = selector.width

                while sel_x.get_value() != x1:
                    reset_say()
                    self.play(self.turn_on(self.grid_coords[(sel_x.get_value(),sel_y.get_value())], fast))
                    D = error_value.get_value()
                    if D>0:
                        say(f"D > 0 \\rightarrow y = y {yis}")
                        self.play(AnimationGroup(
                            error_value.animate.set_value(D + (2*(dy - dx))),
                            sel_y.animate.set_value(sel_y.get_value() + yi),
                            sel_x.animate.set_value(sel_x.get_value() + 1),
                            selector.animate.shift((sels * DOWN * yi) + (sels * RIGHT))
                        ), run_time=speed)
                    else:
                        say("D \\leq 0 \\rightarrow y=y")
                        self.play(AnimationGroup(
                            error_value.animate.set_value(D + (2*dy)),
                            sel_x.animate.set_value(sel_x.get_value() + 1),
                            selector.animate.shift(sels * RIGHT)
                        ), run_time=speed)
                
                self.play(self.turn_on(self.grid_coords[(sel_x.get_value(),sel_y.get_value())], fast))







            if abs(y1-y0) < abs(x1-x0):
                say("dy < dx: low")
                if x0>x1:
                    say("Ax > Bx: swap")
                    self.play(Rotate(ideal_line), run_time=speed)
                    linelow(x1,y1,x0,y0)
                else:
                    say("Ax \\leq Bx: ok")
                    linelow(x0,y0,x1,y1)
            else:
                say("dy \\geq dx: high")
                if y0>y1:
                    say("Ay > By: swap")
                    self.play(Rotate(ideal_line), run_time=speed)
                    linehigh(x1,y1,x0,y0)
                else:
                    say("Ax \\leq Bx: ok")
                    linehigh(x0,y0,x1,y1)



        except Exception as e:
            raise e
        finally:
            self.play(LaggedStartMap(lambda x: FadeOut(x, shift=LEFT), said+displayed))


    def turn_on_grid(self, x, y, fast=False):
        try:
            return self.turn_on(self.grid_coords[(x,y)], fast)
        except KeyError: return Animation(run_time=0)

    def draw_circle(self, x0, y0, radius, fast=False):
        speed = 0.5 if fast else 1

        said = []
        def say(text):
            if fast: return
            tex = MathTex(text)
            if not said:
                tex.move_to(FullScreenRectangle().get_boundary_point(UP+LEFT)).shift_onto_screen()
            else:
                tex.next_to(said[-1], DOWN).align_to(said[-1], LEFT)
            said.append(tex)
            self.play(FadeIn(tex, shift=UP, run_time=0.5))
        
        def reset_say():
            if fast: return
            self.play(LaggedStartMap(lambda x: FadeOut(x, shift=LEFT), said, run_time=0.3))   
            said.clear() 


        center: Square = self.grid_coords[(x0, y0)]
        edge: Square = self.grid_coords.get((x0 + radius, y0)) or self.grid_coords.get((x0 - radius, y0)) or self.grid_coords.get((x0, y0+radius)) or self.grid_coords.get((x0, y0-radius))
        if edge is None: raise Exception("Circle is too big")

        say(f"C = ({x0}, {y0})")
        centerdot = Dot(center.get_center())
        self.play(FadeIn(centerdot, run_time=0.2*speed))
        self.play(Circumscribe(centerdot), run_time=speed)

        say(f"r = {radius}")
        radius_arrow = Line(center.get_center(), edge.get_center())
        self.play(Create(radius_arrow), run_time=speed)

        screen_size = np.linalg.norm(edge.get_center() - center.get_center())

        ideal_circle = Circle(screen_size, WHITE).move_arc_center_to(center.get_center())
        self.play(Create(ideal_circle), run_time=speed)
        self.play(FadeOut(centerdot, radius_arrow, shift=LEFT), run_time=speed)
        displayed = [ideal_circle]


        try:
            sel_x = Integer(radius)
            lx = Tex("x = ")
            sel_y = Integer(0)
            ly = Tex("y = ")
            sel_coord = VGroup(lx, sel_x, ly, sel_y).arrange().shift(3*UP)
            displayed.append(sel_coord)
            self.play(Write(sel_coord), run_time=speed)

            error_label = Tex("Error = ").shift(3*DOWN)
            error_value = Integer(0).next_to(error_label)
            error = VGroup(error_label, error_value)
            displayed.append(error)
            self.play(Write(error), run_time=speed)

            selector: Square = edge.copy().set_fill(RED, 0.6)
            displayed.append(selector)
            self.play(Create(selector), run_time=speed)
            sels = selector.width


            while sel_x.get_value() >= sel_y.get_value():
                reset_say()
                x = sel_x.get_value()
                y = sel_y.get_value()
                anims = [
                    self.turn_on_grid(x0 + x, y0+y, fast),
                    self.turn_on_grid(x0 - x, y0+y, fast),

                    self.turn_on_grid(x0 + y, y0+x, fast),
                    self.turn_on_grid(x0 - y, y0+x, fast),

                    self.turn_on_grid(x0 + x, y0-y, fast),
                    self.turn_on_grid(x0 - x, y0-y, fast),

                    self.turn_on_grid(x0 + y, y0-x, fast),
                    self.turn_on_grid(x0 - y, y0-x, fast),
                ]
                self.play(LaggedStart(*[i for i in anims if i]))

                self.play(AnimationGroup(
                    sel_y.animate.set_value(y + 1),
                    selector.animate.shift(sels * UP),
                    error_value.animate.set_value(error_value.get_value() + 1 + 2*(y+1))
                ), run_time=speed)
                say("2 * (E-x) + 1 > 0?")
                if 2 * (error_value.get_value() - x) + 1 > 0:
                    say("\\rightarrow x = x-1; E = E+1 - 2\\cdot x")
                    self.play(AnimationGroup(
                        sel_x.animate.set_value(x - 1),
                        selector.animate.shift(sels * LEFT),
                        error_value.animate.set_value(error_value.get_value() + (1 - 2*(x-1)))
                    ), run_time=speed)
            return ideal_circle


        except Exception as e:
            raise e
        finally:
            self.play(LaggedStartMap(lambda x: FadeOut(x, shift=LEFT), said+displayed))


    def snap_to_grid(self, x, y):
        p = x*RIGHT + y*UP
        sq = None
        dist = float('inf')
        for square in self.grid:
            nd = np.linalg.norm(square.get_center() - p)
            if nd < dist:
                sq = square
                dist = nd
        return sq


    def draw_line_preview(self, x0, y0, x1, y1, fast=False):
        self.play(Create(Line(
            self.grid_coords[(x0,y0)].get_center(),
            self.grid_coords[(x1,y1)].get_center(),
            color=RED
        )))

    def draw_circle_preview(self, x0, y0, radius, fast=False):
        center = self.grid_coords[(x0,y0)].get_center()
        edge: Square = self.grid_coords.get((x0 + radius, y0)) or self.grid_coords.get((x0 - radius, y0)) or self.grid_coords.get((x0, y0+radius)) or self.grid_coords.get((x0, y0-radius))

        c = Circle(
            np.linalg.norm(edge.get_center() - center)
        ).move_arc_center_to(center)
        self.play(Create(c))
        return c


    def construct(self):
        self.prepare_grid()
        self.wait()

        self.draw_circle(71, 40-1, 18, fast=False)
        c = self.draw_circle(71, 40-1, 24, fast=True)  # Big circle
        #self.draw_line(0,0,142//2,79//2, fast=True)
        
        #print(*self.grid_coords.keys(), sep='\n')
        
        rays = 24
        lines = [
            lambda:self.draw_line(58,40, 84,40), # horizontal
            lambda:self.draw_line(84,40, 71,47), # right front
            lambda:self.draw_line(58,40, 70,24), # left edge
            lambda:self.draw_line(70,47, 70,24), # left tetra line
            lambda:self.draw_line(84,40, 71,24), # right edge
            lambda:self.draw_line(71,47, 71,24), # right tetra line
            lambda:self.draw_line(59,40, 70,47), # left front
        ]


        ctr = 0
        for angle in np.linspace(0, 2*PI, rays):
            p = c.point_at_angle(angle)
            vec = (p - c.get_center())
            vec = vec / np.linalg.norm(vec)
            p = p + (vec * 0.5) 
            next_p = p + (vec*0.5)
            s1 = self.snap_to_grid(p[0], p[1])
            s2 = self.snap_to_grid(next_p[0], next_p[1])

            s1c = s2c = None
            for k,v in self.grid_coords.items():
                if v is s1:
                    s1c = k
                if v is s2:
                    s2c = k

            cs = []
            cs.extend(s1c)
            cs.extend(s2c)
            self.draw_line(*cs, fast=True)

            if ctr%3 == 0:
                try:
                    f = lines.pop(0)
                    f()
                except: pass
            ctr += 1


        while lines:
            try:
                f = lines.pop(0)
                f()
            except: pass

        self.wait(5)
