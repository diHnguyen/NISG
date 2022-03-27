global arcs = [1 21; 1 29; 1 30; 1 34; 1 35; 2 3; 2 9; 2 10; 2 14; 2 28; 2 30; 3 4; 3 8; 3 10; 3 22; 3 26; 4 5; 4 32; 4 39; 5 12; 5 14; 5 25; 6 3; 6 29; 6 31; 6 36; 6 38; 7 8; 7 12; 7 19; 7 25; 8 13; 8 15; 8 20; 8 21; 8 24; 9 7; 9 13; 9 20; 9 23; 9 39; 9 40; 10 7; 10 19; 10 28; 10 30; 10 37; 11 8; 11 9; 11 13; 11 21; 11 24; 11 25; 11 31; 11 39; 12 22; 12 25; 12 27; 13 4; 13 7; 13 21; 13 27; 14 6; 14 17; 14 33; 15 16; 15 27; 16 2; 16 3; 16 9; 16 13; 16 27; 16 29; 16 36; 17 20; 17 34; 17 37; 18 2; 18 8; 18 9; 18 25; 18 38; 19 4; 19 14; 19 25; 19 31; 20 4; 20 26; 20 35; 21 5; 21 7; 21 14; 21 19; 21 23; 21 25; 22 2; 22 3; 22 6; 22 12; 23 25; 24 7; 24 16; 24 20; 24 27; 25 10; 25 18; 25 40; 26 12; 26 17; 26 21; 26 25; 26 29; 27 4; 27 15; 27 17; 27 31; 28 10; 28 12; 28 20; 28 25; 28 26; 28 30; 28 40; 29 17; 29 25; 30 10; 30 15; 30 17; 30 28; 30 31; 30 39; 31 2; 31 4; 31 14; 31 25; 31 28; 32 22; 33 3; 33 11; 33 19; 33 25; 34 11; 34 18; 34 20; 34 26; 34 31; 35 13; 35 23; 35 38; 36 13; 36 32; 36 38; 36 39; 36 40; 37 2; 37 4; 37 32; 37 39; 38 33; 38 35; 39 4; 39 9; 39 14; 39 16; 39 22]
global d_x = [6.0, 2.0, 6.0, 3.0, 6.0, 4.0, 4.0, 8.0, 10.0, 1.0, 9.0, 9.0, 10.0, 3.0, 6.0, 2.0, 10.0, 2.0, 1.0, 8.0, 10.0, 5.0, 2.0, 3.0, 6.0, 7.0, 5.0, 7.0, 4.0, 7.0, 9.0, 4.0, 9.0, 8.0, 3.0, 4.0, 2.0, 9.0, 9.0, 8.0, 3.0, 8.0, 7.0, 1.0, 2.0, 7.0, 4.0, 3.0, 8.0, 5.0, 3.0, 9.0, 4.0, 2.0, 10.0, 4.0, 6.0, 10.0, 5.0, 4.0, 9.0, 5.0, 6.0, 6.0, 2.0, 4.0, 5.0, 7.0, 4.0, 6.0, 8.0, 4.0, 10.0, 2.0, 2.0, 5.0, 9.0, 9.0, 3.0, 6.0, 7.0, 1.0, 3.0, 6.0, 9.0, 6.0, 7.0, 8.0, 1.0, 7.0, 8.0, 9.0, 5.0, 9.0, 2.0, 2.0, 4.0, 4.0, 9.0, 4.0, 4.0, 9.0, 5.0, 10.0, 3.0, 5.0, 5.0, 4.0, 2.0, 2.0, 3.0, 1.0, 4.0, 8.0, 9.0, 5.0, 9.0, 5.0, 9.0, 10.0, 8.0, 2.0, 8.0, 3.0, 8.0, 3.0, 9.0, 7.0, 2.0, 2.0, 6.0, 3.0, 8.0, 1.0, 1.0, 10.0, 2.0, 7.0, 8.0, 5.0, 2.0, 1.0, 10.0, 1.0, 9.0, 1.0, 6.0, 5.0, 3.0, 2.0, 3.0, 4.0, 5.0, 10.0, 10.0, 6.0, 10.0, 4.0, 7.0, 4.0, 8.0, 9.0, 3.0, 3.0, 9.0]
global b_x = 5
global d_y = [7.0, 2.0, 5.0, 10.0, 6.0, 7.0, 10.0, 2.0, 4.0, 9.0, 6.0, 2.0, 5.0, 1.0, 2.0, 6.0, 8.0, 8.0, 2.0, 10.0, 9.0, 8.0, 3.0, 7.0, 6.0, 4.0, 4.0, 1.0, 2.0, 6.0, 1.0, 5.0, 10.0, 4.0, 6.0, 3.0, 6.0, 1.0, 8.0, 3.0, 4.0, 6.0, 3.0, 3.0, 8.0, 6.0, 9.0, 1.0, 6.0, 4.0, 1.0, 5.0, 3.0, 6.0, 5.0, 7.0, 3.0, 10.0, 6.0, 8.0, 5.0, 2.0, 9.0, 1.0, 1.0, 6.0, 2.0, 7.0, 3.0, 7.0, 8.0, 9.0, 9.0, 10.0, 7.0, 10.0, 6.0, 6.0, 8.0, 8.0, 1.0, 9.0, 7.0, 2.0, 4.0, 5.0, 10.0, 9.0, 2.0, 6.0, 8.0, 10.0, 9.0, 8.0, 10.0, 1.0, 4.0, 4.0, 4.0, 10.0, 4.0, 5.0, 9.0, 1.0, 6.0, 7.0, 6.0, 9.0, 6.0, 9.0, 4.0, 4.0, 10.0, 4.0, 10.0, 5.0, 8.0, 3.0, 3.0, 1.0, 1.0, 5.0, 3.0, 10.0, 7.0, 8.0, 6.0, 5.0, 7.0, 5.0, 1.0, 1.0, 6.0, 5.0, 10.0, 8.0, 8.0, 4.0, 1.0, 9.0, 2.0, 4.0, 4.0, 1.0, 2.0, 2.0, 1.0, 9.0, 8.0, 6.0, 2.0, 3.0, 7.0, 10.0, 4.0, 9.0, 8.0, 6.0, 9.0, 10.0, 10.0, 10.0, 9.0, 4.0, 6.0]
global b_y = 10
global p = [0.112, 0.004, 0.216, 0.516, 0.955, 0.568, 0.285, 0.225, 0.73, 0.922, 0.255, 0.47, 0.993, 0.593, 0.221, 0.617, 0.382, 0.36, 0.788, 0.721, 0.224, 0.246, 0.735, 0.543, 0.737, 0.14, 0.707, 0.19, 0.842, 0.72, 0.78, 0.956, 0.553, 0.184, 0.03, 0.656, 0.668, 0.34, 0.715, 0.767, 0.197, 0.133, 0.703, 0.261, 0.243, 0.97, 0.717, 0.947, 0.212, 0.995, 0.921, 0.172, 0.26, 0.288, 0.215, 0.406, 0.662, 0.285, 0.764, 0.85, 0.846, 0.032, 0.433, 0.568, 0.636, 0.685, 0.564, 0.153, 0.536, 0.573, 0.999, 0.9, 0.92, 0.455, 0.285, 0.896, 0.872, 0.966, 0.046, 0.925, 0.038, 0.718, 0.432, 0.203, 0.796, 0.429, 0.188, 0.209, 0.182, 0.316, 0.832, 0.7, 0.907, 0.073, 0.649, 0.58, 0.988, 0.569, 0.08, 0.095, 0.221, 0.105, 0.564, 0.468, 0.964, 0.925, 0.383, 0.127, 0.366, 0.089, 0.639, 0.425, 0.673, 0.984, 0.955, 0.276, 0.008, 0.35, 0.362, 0.356, 0.709, 0.978, 0.218, 0.562, 0.232, 0.273, 0.593, 0.665, 0.758, 0.495, 0.4, 0.787, 0.763, 0.983, 0.048, 0.566, 0.075, 0.46, 0.784, 0.968, 0.818, 0.051, 0.154, 0.279, 0.056, 0.532, 0.939, 0.429, 0.429, 0.378, 0.716, 0.276, 0.352, 0.309, 0.792, 0.683, 0.161, 0.385, 0.195, 0.671, 0.567, 0.053, 0.447, 0.373, 0.258]
global q = [0.209, 0.962, 0.531, 0.531, 0.988, 0.857, 0.509, 0.503, 0.914, 0.959, 0.958, 0.663, 0.997, 0.869, 0.51, 0.886, 0.991, 0.739, 0.888, 0.934, 0.755, 0.394, 0.975, 0.597, 0.744, 0.586, 0.931, 0.841, 0.931, 0.873, 0.853, 0.978, 0.85, 0.341, 0.803, 0.672, 0.695, 0.763, 0.771, 0.859, 0.297, 0.236, 0.949, 0.835, 0.63, 0.989, 0.725, 0.987, 0.879, 0.996, 0.938, 0.175, 0.549, 0.342, 0.546, 0.724, 0.837, 0.877, 0.914, 0.94, 0.995, 0.331, 0.654, 0.704, 0.99, 0.735, 0.668, 0.651, 0.679, 0.679, 0.999, 0.977, 0.966, 0.911, 0.756, 0.998, 0.982, 0.995, 0.706, 0.965, 0.585, 0.85, 0.683, 0.814, 0.806, 0.988, 0.973, 0.91, 0.291, 0.737, 0.878, 0.954, 0.982, 0.983, 0.797, 0.894, 0.999, 0.973, 0.316, 0.636, 0.867, 0.99, 0.643, 0.532, 0.978, 0.995, 0.911, 0.303, 0.922, 0.158, 0.939, 0.652, 0.974, 0.989, 0.957, 0.794, 0.577, 0.953, 0.466, 0.465, 0.871, 0.997, 0.833, 0.908, 0.712, 0.808, 0.799, 0.763, 0.92, 0.881, 0.492, 0.96, 0.955, 0.989, 0.231, 0.846, 0.378, 0.584, 0.814, 0.98, 0.891, 0.452, 0.52, 0.807, 0.974, 0.9, 0.957, 0.952, 0.705, 0.508, 0.906, 0.566, 0.7, 0.436, 0.986, 0.978, 0.232, 0.577, 0.788, 0.757, 0.863, 0.165, 0.667, 0.991, 0.266]
global origin = 1
global destination = 40