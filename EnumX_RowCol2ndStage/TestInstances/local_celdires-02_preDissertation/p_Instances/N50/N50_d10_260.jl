global arcs = [1 10; 1 27; 1 34; 1 35; 1 45; 1 46; 1 48; 2 14; 2 40; 3 2; 3 20; 3 44; 3 48; 4 9; 4 20; 4 34; 4 45; 5 24; 5 40; 5 49; 5 50; 6 14; 6 16; 6 26; 6 40; 6 47; 7 13; 7 29; 7 32; 7 34; 7 40; 7 47; 8 6; 8 29; 8 32; 8 41; 8 45; 9 6; 9 12; 9 18; 9 30; 9 40; 9 47; 10 11; 10 16; 10 37; 10 41; 10 47; 11 4; 11 17; 11 26; 11 27; 11 30; 11 36; 11 40; 12 13; 12 15; 12 26; 12 37; 13 14; 13 21; 13 26; 13 28; 13 33; 13 47; 14 3; 14 6; 15 11; 15 17; 15 44; 15 45; 16 11; 16 20; 16 26; 16 27; 16 45; 17 4; 17 7; 17 35; 17 37; 17 45; 18 7; 18 24; 19 18; 19 33; 19 38; 19 43; 20 6; 20 11; 20 19; 20 39; 21 11; 22 15; 22 33; 22 44; 22 50; 23 17; 24 4; 24 9; 24 13; 24 20; 24 30; 24 31; 24 35; 24 40; 24 46; 25 3; 25 4; 25 10; 25 14; 25 32; 25 45; 26 12; 26 18; 26 19; 26 22; 26 38; 26 40; 27 8; 27 16; 27 17; 27 37; 27 46; 27 49; 28 30; 28 33; 28 41; 29 22; 29 37; 30 4; 30 12; 30 23; 30 41; 30 45; 31 7; 31 19; 31 24; 31 30; 31 42; 31 50; 32 3; 32 12; 32 28; 32 35; 32 50; 33 3; 33 10; 33 19; 33 20; 33 28; 33 34; 33 43; 33 45; 34 2; 34 13; 34 19; 34 20; 34 32; 35 9; 35 14; 35 40; 36 3; 36 4; 36 20; 36 42; 37 4; 37 12; 38 4; 38 5; 38 10; 38 13; 38 16; 38 31; 38 36; 38 40; 39 33; 39 38; 39 40; 39 49; 40 12; 40 14; 40 18; 40 23; 40 26; 41 24; 41 25; 41 36; 41 46; 42 9; 42 17; 42 26; 42 40; 42 45; 43 16; 43 25; 43 31; 43 38; 43 41; 44 7; 44 25; 44 30; 44 49; 45 14; 45 15; 45 20; 45 41; 45 48; 46 12; 46 29; 46 31; 46 35; 46 40; 46 45; 46 50; 47 8; 47 11; 47 12; 47 42; 48 18; 48 19; 48 33; 48 34; 48 43; 48 44; 49 7; 49 15; 49 17; 49 36; 49 38]
global d_x = [10.0, 3.0, 9.0, 7.0, 7.0, 5.0, 10.0, 7.0, 8.0, 4.0, 2.0, 10.0, 3.0, 1.0, 2.0, 7.0, 3.0, 7.0, 6.0, 8.0, 7.0, 1.0, 3.0, 7.0, 6.0, 8.0, 9.0, 2.0, 4.0, 8.0, 10.0, 9.0, 7.0, 6.0, 9.0, 6.0, 4.0, 6.0, 9.0, 9.0, 3.0, 2.0, 1.0, 9.0, 8.0, 9.0, 2.0, 10.0, 10.0, 5.0, 10.0, 5.0, 1.0, 1.0, 4.0, 4.0, 6.0, 1.0, 1.0, 9.0, 3.0, 7.0, 8.0, 4.0, 5.0, 3.0, 9.0, 7.0, 3.0, 3.0, 4.0, 6.0, 8.0, 10.0, 10.0, 7.0, 8.0, 4.0, 3.0, 3.0, 9.0, 10.0, 10.0, 4.0, 5.0, 3.0, 7.0, 2.0, 6.0, 1.0, 5.0, 7.0, 7.0, 6.0, 10.0, 1.0, 7.0, 4.0, 2.0, 7.0, 8.0, 3.0, 9.0, 7.0, 3.0, 2.0, 6.0, 1.0, 1.0, 2.0, 4.0, 4.0, 6.0, 6.0, 1.0, 10.0, 8.0, 3.0, 8.0, 1.0, 4.0, 4.0, 5.0, 6.0, 6.0, 9.0, 6.0, 1.0, 3.0, 1.0, 10.0, 6.0, 3.0, 9.0, 4.0, 10.0, 9.0, 5.0, 9.0, 1.0, 7.0, 2.0, 1.0, 7.0, 3.0, 7.0, 3.0, 9.0, 4.0, 8.0, 8.0, 8.0, 10.0, 7.0, 7.0, 2.0, 3.0, 8.0, 9.0, 4.0, 1.0, 5.0, 9.0, 1.0, 3.0, 2.0, 10.0, 10.0, 10.0, 2.0, 7.0, 10.0, 2.0, 10.0, 9.0, 7.0, 10.0, 5.0, 5.0, 3.0, 2.0, 3.0, 10.0, 2.0, 1.0, 8.0, 4.0, 7.0, 8.0, 7.0, 4.0, 5.0, 8.0, 3.0, 4.0, 3.0, 6.0, 6.0, 5.0, 8.0, 6.0, 1.0, 8.0, 1.0, 4.0, 8.0, 6.0, 7.0, 3.0, 5.0, 10.0, 9.0, 8.0, 10.0, 4.0, 3.0, 2.0, 9.0, 6.0, 10.0, 5.0, 9.0, 7.0, 2.0, 7.0, 2.0, 10.0, 9.0, 1.0]
global b_x = 5
global d_y = [8.0, 3.0, 9.0, 9.0, 10.0, 4.0, 10.0, 4.0, 5.0, 6.0, 6.0, 9.0, 8.0, 7.0, 5.0, 7.0, 9.0, 7.0, 3.0, 1.0, 3.0, 3.0, 7.0, 10.0, 9.0, 5.0, 3.0, 6.0, 1.0, 3.0, 2.0, 7.0, 4.0, 7.0, 9.0, 8.0, 5.0, 2.0, 1.0, 8.0, 1.0, 5.0, 4.0, 4.0, 6.0, 5.0, 3.0, 3.0, 7.0, 3.0, 3.0, 5.0, 4.0, 6.0, 6.0, 8.0, 1.0, 4.0, 6.0, 7.0, 1.0, 1.0, 9.0, 7.0, 10.0, 8.0, 6.0, 10.0, 1.0, 3.0, 10.0, 8.0, 1.0, 9.0, 1.0, 4.0, 8.0, 2.0, 9.0, 5.0, 6.0, 3.0, 2.0, 4.0, 6.0, 10.0, 2.0, 8.0, 3.0, 4.0, 3.0, 2.0, 1.0, 9.0, 8.0, 3.0, 2.0, 9.0, 5.0, 10.0, 3.0, 1.0, 2.0, 5.0, 4.0, 8.0, 7.0, 5.0, 1.0, 1.0, 9.0, 5.0, 1.0, 8.0, 7.0, 9.0, 1.0, 5.0, 6.0, 5.0, 10.0, 3.0, 7.0, 5.0, 8.0, 2.0, 9.0, 9.0, 5.0, 4.0, 2.0, 4.0, 9.0, 8.0, 4.0, 3.0, 1.0, 2.0, 6.0, 6.0, 9.0, 2.0, 7.0, 10.0, 4.0, 4.0, 6.0, 8.0, 6.0, 4.0, 9.0, 5.0, 4.0, 2.0, 5.0, 2.0, 10.0, 3.0, 7.0, 5.0, 7.0, 4.0, 9.0, 2.0, 4.0, 1.0, 10.0, 2.0, 2.0, 1.0, 8.0, 9.0, 4.0, 8.0, 2.0, 1.0, 5.0, 7.0, 5.0, 3.0, 6.0, 6.0, 10.0, 10.0, 6.0, 7.0, 5.0, 7.0, 5.0, 9.0, 8.0, 7.0, 7.0, 9.0, 1.0, 5.0, 6.0, 7.0, 8.0, 5.0, 9.0, 5.0, 9.0, 8.0, 7.0, 2.0, 7.0, 5.0, 4.0, 3.0, 4.0, 6.0, 6.0, 8.0, 6.0, 9.0, 8.0, 8.0, 9.0, 1.0, 9.0, 7.0, 3.0, 7.0, 2.0, 1.0, 8.0, 5.0, 4.0]
global b_y = 10
global p = [0.708, 0.76, 0.454, 0.933, 0.661, 0.789, 0.867, 0.1, 0.192, 0.417, 0.179, 0.571, 0.863, 0.244, 0.575, 0.372, 0.396, 0.411, 0.462, 0.413, 0.48, 0.836, 0.509, 0.428, 0.958, 0.666, 0.314, 0.534, 0.187, 0.239, 0.308, 0.012, 0.926, 0.818, 0.134, 0.745, 0.818, 0.295, 0.82, 0.808, 0.939, 0.309, 0.264, 0.775, 0.063, 0.113, 0.837, 0.872, 0.424, 0.555, 0.482, 0.128, 0.345, 0.099, 0.512, 0.243, 0.681, 0.626, 0.863, 0.853, 0.064, 0.59, 0.358, 0.496, 0.762, 0.764, 0.797, 0.266, 0.499, 0.264, 0.767, 0.307, 0.356, 0.171, 0.377, 0.487, 0.721, 0.772, 0.666, 0.213, 0.527, 0.328, 0.398, 0.444, 0.737, 0.957, 0.911, 0.702, 0.992, 0.919, 0.744, 0.769, 0.813, 0.165, 0.62, 0.864, 0.871, 0.384, 0.67, 0.909, 0.655, 0.748, 0.316, 0.219, 0.497, 0.471, 0.334, 0.279, 0.395, 0.301, 0.902, 0.383, 0.399, 0.073, 0.748, 0.85, 0.393, 0.549, 0.616, 0.247, 0.248, 0.16, 0.281, 0.813, 0.512, 0.053, 0.06, 0.793, 0.686, 0.235, 0.748, 0.074, 0.029, 0.002, 0.365, 0.802, 0.749, 0.162, 0.17, 0.849, 0.377, 0.382, 0.628, 0.081, 0.593, 0.59, 0.71, 0.813, 0.452, 0.243, 0.101, 0.859, 0.953, 0.971, 0.456, 0.225, 0.736, 0.085, 0.52, 0.238, 0.871, 0.819, 0.211, 0.864, 0.087, 0.284, 0.205, 0.804, 0.523, 0.785, 0.467, 0.67, 0.294, 0.012, 0.59, 0.693, 0.073, 0.174, 0.572, 0.619, 0.903, 0.477, 0.06, 0.509, 0.034, 0.886, 0.123, 0.767, 0.606, 0.537, 0.264, 0.045, 0.399, 0.695, 0.736, 0.808, 0.121, 0.742, 0.742, 0.641, 0.877, 0.45, 0.665, 0.371, 0.623, 0.45, 0.681, 0.762, 0.976, 0.143, 0.16, 0.205, 0.605, 0.922, 0.641, 0.837, 0.444, 0.779, 0.609, 0.148, 0.042, 0.167, 0.482, 0.094, 0.569, 0.124, 0.056, 0.992, 0.675]
global q = [0.865, 0.768, 0.756, 0.987, 0.995, 0.898, 0.915, 0.694, 0.84, 0.633, 0.922, 0.778, 0.944, 0.909, 0.997, 0.952, 0.596, 0.477, 0.573, 0.926, 0.83, 0.926, 0.562, 0.69, 0.993, 0.866, 0.67, 0.817, 0.245, 0.462, 0.518, 0.067, 0.987, 0.977, 0.684, 0.775, 0.955, 0.909, 0.901, 0.878, 0.951, 0.854, 0.498, 0.829, 0.281, 0.596, 0.962, 0.963, 0.67, 0.975, 0.638, 0.505, 0.95, 0.817, 0.642, 0.467, 0.844, 0.776, 0.885, 0.981, 0.349, 0.82, 0.685, 0.713, 0.931, 0.812, 0.803, 0.645, 0.988, 0.416, 0.892, 0.37, 0.902, 0.74, 0.404, 0.704, 0.926, 0.943, 0.687, 0.46, 0.79, 0.759, 0.862, 0.883, 0.942, 0.98, 0.957, 0.756, 0.993, 0.972, 0.904, 0.95, 0.968, 0.591, 0.793, 0.944, 0.889, 0.553, 0.945, 0.984, 0.698, 0.917, 0.996, 0.66, 0.894, 0.873, 0.791, 0.716, 0.927, 0.567, 0.979, 0.675, 0.996, 0.083, 0.975, 0.878, 0.976, 0.589, 0.954, 0.982, 0.904, 0.682, 0.56, 0.971, 0.86, 0.166, 0.878, 0.899, 0.721, 0.968, 0.988, 0.465, 0.807, 0.944, 0.472, 0.83, 0.954, 0.62, 0.959, 0.886, 0.727, 0.798, 0.81, 0.626, 0.656, 0.677, 0.842, 0.916, 0.586, 0.731, 0.822, 0.987, 0.988, 0.985, 0.488, 0.546, 0.784, 0.909, 0.826, 0.9, 0.966, 0.959, 0.425, 0.883, 0.222, 0.376, 0.523, 0.86, 0.606, 0.873, 0.549, 0.917, 0.882, 0.09, 0.648, 0.889, 0.533, 0.333, 0.811, 0.619, 0.981, 0.757, 0.178, 0.693, 0.76, 0.91, 0.691, 0.828, 0.627, 0.865, 0.507, 0.594, 0.625, 0.733, 0.844, 0.922, 0.798, 0.99, 0.761, 0.801, 0.908, 0.606, 0.951, 0.739, 0.787, 0.596, 0.766, 0.879, 0.995, 0.309, 0.361, 0.403, 0.854, 0.942, 0.824, 0.956, 0.862, 0.887, 0.736, 0.291, 0.827, 0.887, 0.804, 0.874, 0.595, 0.823, 0.111, 0.993, 0.85]
global origin = 1
global destination = 50