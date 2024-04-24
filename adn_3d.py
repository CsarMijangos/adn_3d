from manim import *
import numpy as np

class ADN_3D(ThreeDScene):
    def construct(self):
        
        axis_config = {
            "x_range": (-8,8,1),
            "y_range": (-8,8,1),
            "z_range": (-8,8,1)
        }
        nucleotides_config = {
            "radio" : 0.2,
            "color_A": BLUE,
            "color_C": GREEN,
            "color_G": YELLOW,
            "color_T": RED    
        }

        nucle_colors = {
            "A": BLUE,
            "C": GREEN,
            "G" : YELLOW,
            "T" : RED
        }

        helix_config = {
            "length_range": [-2*PI, 2*PI], # extremo derecho solo multiplos de 2*PI
            "color": YELLOW
        }

        puente_config = {
            "thickness" : 0.05,
            "color": WHITE
        }

        complements = {
            "A":"T",
            "C": "G",
            "G": "C",
            "T": "A"
        }


        sequence = ["A", "A", "T", "G", "C", "T", "A", "G", "C", "G", "A", "A", "C", "G", "C", "T", "A", "T", "T", "C"]
        #axes = ThreeDAxes(**axis_config)
        self.set_camera_orientation(phi=80*DEGREES, theta=-40*DEGREES, distance=6)

        #text3d = Text("Hello 3D world!").scale(2)
        #text3d.rotate(PI/2, axis= RIGHT)

        #self.play(Write(axes))

        #self.play(Write(text3d))

        ## Dibujar hebras:
        
        #curve_1 = ParametricFunction(lambda t: np.array([np.cos(t), np.sin(t), t]), 
        #                          color = helix_config["color"], t_range= helix_config["length_range"])
        
        #curve_2 = ParametricFunction(lambda t: np.array([-np.cos(t), -np.sin(t), t]), 
        #                          color = helix_config["color"], t_range= helix_config["length_range"])
        
        def nucleotide_position(seq_list):
            """ This function returns the position of the nucleotides over the helix_curve. 
            It fitts 10 nucleotides by twist.

            Args:
                seq_list (list): A list of A, C, G or T to draw in curve 1
                helix_len (integer): Multiples of 2*PI that represent the length of the helix
            """

            complement_seq = [complements[x] for x in seq_list]

            hebra_1 = VGroup()
            hebra_2 = VGroup()
            puentes = VGroup()

            delta = 9 # uno menos que el número de pb por periodo 
            for i in range(0,len(seq_list)):
                hebra_1.add(Dot3D(point=np.array([0.5*np.cos(-2*PI + 2*PI*i/delta), 0.5*np.sin(-2*PI + 2*PI*i/delta),-2*PI + 2*PI*i/delta]),
                                  radius=nucleotides_config["radio"], 
                                  color=nucle_colors[seq_list[i]],
                                  resolution=(12, 12)))
                hebra_2.add(Dot3D(point=np.array([-0.5*np.cos(-2*PI + 2*PI*i/delta), -0.5*np.sin(-2*PI + 2*PI*i/delta), -2*PI + 2*PI*i/delta]),
                                  radius=nucleotides_config["radio"], 
                                  color=nucle_colors[complement_seq[i]],
                                  resolution=(12, 12)))
                
                puentes.add(Line3D(start=hebra_1[i],end= hebra_2[i],
                                   thickness=puente_config["thickness"], 
                                   color=puente_config["color"]))
            
            #return VGroup(hebra_1, hebra_2)
            return VGroup(hebra_1, hebra_2, puentes)


        #nuc_a = Dot3D(point=np.array([np.cos(0), np.sin(0), 0]),radius=nucleotides_config["radio"], color=nucleotides_config["color_A"])
        #nuc_t = Dot3D(point=np.array([-np.cos(0), -np.sin(0), -0]),radius=nucleotides_config["radio"], color=nucleotides_config["color_T"])
        #bridge = Line3D(start=nuc_a,end= nuc_t, thickness=0.02, color=WHITE)

        #step = VGroup(nuc_a,nuc_t, bridge)

        #self.play(Write(step))

        #hebras = VGroup(curve_1,curve_2)

        nucleotidos = nucleotide_position(sequence)

        self.begin_ambient_camera_rotation(rate = 1)

        self.play(Write(nucleotidos), run_time=8)
        
        #self.play(Write(hebras), run_time=8)

        #self.add(hebras)
        #self.add(nucleotidos)
        


        self.wait(4)