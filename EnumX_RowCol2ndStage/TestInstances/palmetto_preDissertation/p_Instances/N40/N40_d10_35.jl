global arcs = [1 11; 1 16; 1 30; 1 34; 1 37; 2 22; 2 30; 2 33; 3 12; 3 18; 3 30; 4 17; 4 39; 5 2; 5 8; 5 15; 5 24; 6 15; 6 26; 7 6; 7 20; 7 30; 8 12; 8 34; 9 5; 9 20; 9 24; 9 34; 9 37; 9 39; 10 6; 10 19; 11 5; 11 9; 11 27; 11 33; 11 35; 11 40; 12 23; 12 24; 12 32; 12 40; 13 2; 13 8; 13 25; 13 26; 13 36; 14 13; 14 16; 14 19; 14 20; 14 24; 14 32; 15 12; 15 16; 15 23; 15 28; 15 35; 15 39; 16 11; 17 5; 17 11; 17 32; 17 39; 18 5; 18 11; 18 28; 19 5; 20 10; 20 17; 20 18; 20 26; 20 38; 21 7; 21 16; 21 39; 22 3; 22 4; 22 12; 22 18; 22 29; 22 31; 22 33; 23 7; 23 24; 23 39; 24 15; 25 16; 26 14; 26 34; 26 38; 27 4; 27 12; 27 19; 27 22; 27 28; 28 3; 28 18; 28 29; 28 37; 29 22; 29 23; 30 2; 30 19; 30 20; 30 31; 31 13; 31 20; 31 39; 32 10; 32 29; 33 3; 33 10; 33 15; 33 21; 33 24; 33 40; 34 5; 34 18; 35 3; 35 25; 36 11; 36 13; 37 11; 37 14; 37 36; 38 19; 38 20; 38 21; 38 31; 38 33; 38 40; 39 8; 39 9; 39 19]
global d_x = [5.0, 2.0, 10.0, 3.0, 1.0, 1.0, 9.0, 7.0, 10.0, 10.0, 2.0, 1.0, 1.0, 6.0, 1.0, 8.0, 9.0, 7.0, 1.0, 3.0, 4.0, 7.0, 4.0, 4.0, 7.0, 2.0, 7.0, 9.0, 2.0, 7.0, 1.0, 2.0, 1.0, 4.0, 9.0, 3.0, 2.0, 7.0, 4.0, 2.0, 1.0, 8.0, 5.0, 4.0, 2.0, 6.0, 2.0, 3.0, 7.0, 1.0, 4.0, 10.0, 10.0, 1.0, 1.0, 1.0, 8.0, 7.0, 4.0, 4.0, 3.0, 6.0, 10.0, 5.0, 8.0, 10.0, 9.0, 5.0, 6.0, 2.0, 1.0, 1.0, 2.0, 1.0, 4.0, 8.0, 6.0, 3.0, 6.0, 6.0, 4.0, 5.0, 8.0, 2.0, 6.0, 2.0, 10.0, 1.0, 9.0, 6.0, 3.0, 10.0, 5.0, 8.0, 7.0, 2.0, 5.0, 6.0, 4.0, 10.0, 10.0, 3.0, 3.0, 7.0, 6.0, 2.0, 8.0, 8.0, 10.0, 7.0, 8.0, 9.0, 8.0, 10.0, 3.0, 10.0, 7.0, 3.0, 8.0, 8.0, 2.0, 8.0, 8.0, 8.0, 3.0, 3.0, 6.0, 9.0, 10.0, 9.0, 9.0, 1.0, 1.0, 5.0, 5.0]
global b_x = 5
global d_y = [1.0, 3.0, 2.0, 9.0, 7.0, 6.0, 4.0, 7.0, 4.0, 4.0, 4.0, 7.0, 6.0, 7.0, 3.0, 8.0, 7.0, 2.0, 6.0, 9.0, 6.0, 4.0, 6.0, 8.0, 7.0, 1.0, 10.0, 3.0, 5.0, 5.0, 6.0, 7.0, 10.0, 10.0, 1.0, 8.0, 8.0, 8.0, 3.0, 5.0, 6.0, 3.0, 4.0, 7.0, 6.0, 1.0, 2.0, 8.0, 1.0, 5.0, 9.0, 5.0, 5.0, 8.0, 7.0, 4.0, 9.0, 10.0, 5.0, 8.0, 10.0, 4.0, 6.0, 6.0, 2.0, 9.0, 4.0, 9.0, 7.0, 6.0, 2.0, 2.0, 1.0, 7.0, 7.0, 10.0, 4.0, 4.0, 5.0, 8.0, 3.0, 8.0, 5.0, 7.0, 1.0, 5.0, 1.0, 2.0, 5.0, 6.0, 7.0, 7.0, 8.0, 7.0, 8.0, 9.0, 6.0, 9.0, 10.0, 1.0, 2.0, 1.0, 4.0, 5.0, 4.0, 5.0, 9.0, 7.0, 2.0, 6.0, 3.0, 5.0, 9.0, 2.0, 3.0, 8.0, 8.0, 3.0, 8.0, 1.0, 10.0, 4.0, 2.0, 8.0, 1.0, 6.0, 8.0, 8.0, 1.0, 6.0, 10.0, 5.0, 1.0, 2.0, 3.0]
global b_y = 10
global p = [0.349, 0.392, 0.516, 0.282, 0.396, 0.987, 0.366, 0.485, 0.393, 0.679, 0.624, 0.411, 0.613, 0.425, 0.843, 0.296, 0.511, 0.169, 0.632, 0.042, 0.648, 0.378, 0.38, 0.203, 0.645, 0.442, 0.846, 0.928, 0.494, 0.558, 0.812, 0.549, 0.125, 0.472, 0.216, 0.448, 0.085, 0.668, 0.985, 0.363, 0.054, 0.55, 0.004, 0.17, 0.761, 0.176, 0.575, 0.836, 0.951, 0.818, 0.387, 0.327, 0.179, 0.505, 0.241, 0.818, 0.921, 0.779, 0.625, 0.848, 0.913, 0.345, 0.659, 0.232, 0.435, 0.763, 0.943, 0.025, 0.406, 0.619, 0.426, 0.649, 0.484, 0.172, 0.982, 0.523, 0.114, 0.967, 0.371, 0.012, 0.06, 0.177, 0.693, 0.356, 0.629, 0.549, 0.993, 0.247, 0.255, 0.455, 0.897, 0.196, 0.836, 0.569, 0.377, 0.511, 0.631, 0.798, 0.185, 0.52, 0.934, 0.946, 0.647, 0.029, 0.084, 0.475, 0.647, 0.73, 0.538, 0.774, 0.647, 0.557, 0.362, 0.543, 0.294, 0.543, 0.648, 0.471, 0.356, 0.984, 0.486, 0.678, 0.017, 0.453, 0.253, 0.647, 0.783, 0.359, 0.07, 0.43, 0.929, 0.349, 0.972, 0.93, 0.846]
global q = [0.471, 0.797, 0.718, 0.405, 0.42, 0.993, 0.992, 0.668, 0.826, 0.733, 0.645, 0.777, 0.855, 0.616, 0.99, 0.984, 0.966, 0.501, 0.724, 0.464, 0.925, 0.473, 0.669, 0.855, 0.673, 0.819, 0.961, 0.947, 0.844, 0.678, 0.966, 0.739, 0.445, 0.953, 0.776, 0.889, 0.783, 0.93, 0.988, 0.506, 0.09, 0.554, 0.178, 0.307, 0.987, 0.97, 0.731, 0.973, 0.978, 0.972, 0.588, 0.661, 0.632, 0.904, 0.922, 0.858, 0.926, 0.927, 0.908, 0.901, 0.93, 0.577, 0.974, 0.558, 0.952, 0.91, 0.971, 0.693, 0.486, 0.876, 0.994, 0.663, 0.864, 0.744, 0.994, 0.913, 0.279, 0.998, 0.828, 0.645, 0.698, 0.777, 0.949, 0.367, 0.75, 0.608, 0.996, 0.413, 0.695, 0.648, 0.946, 0.234, 0.929, 0.937, 0.537, 0.826, 0.881, 0.982, 0.627, 0.837, 0.974, 0.992, 0.803, 0.271, 0.577, 0.686, 0.731, 0.878, 0.684, 0.888, 0.926, 0.784, 0.542, 0.621, 0.294, 0.55, 0.895, 0.818, 0.558, 0.995, 0.801, 0.99, 0.964, 0.518, 0.385, 0.847, 0.863, 0.704, 0.44, 0.952, 0.98, 0.449, 0.988, 0.992, 0.88]
global origin = 1
global destination = 40