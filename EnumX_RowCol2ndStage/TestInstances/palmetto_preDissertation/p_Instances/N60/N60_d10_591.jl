global arcs = [1 3; 1 10; 1 15; 1 25; 1 27; 1 36; 1 40; 1 44; 1 51; 1 60; 2 13; 2 26; 2 51; 3 2; 3 51; 3 52; 3 54; 4 22; 4 34; 4 35; 4 50; 5 4; 5 6; 5 19; 5 31; 5 33; 5 37; 5 54; 5 55; 6 2; 6 10; 6 11; 6 23; 6 43; 6 46; 6 52; 6 57; 7 16; 7 19; 7 33; 7 34; 7 50; 8 9; 8 41; 8 52; 8 58; 8 60; 9 16; 9 18; 9 42; 9 43; 10 27; 10 31; 10 36; 10 41; 10 55; 11 22; 11 30; 11 44; 11 48; 11 51; 11 59; 12 4; 12 5; 12 9; 12 22; 12 50; 12 56; 13 3; 13 6; 13 17; 13 23; 13 26; 13 32; 13 52; 14 16; 14 22; 14 24; 14 30; 14 33; 14 53; 15 12; 15 44; 15 46; 15 50; 15 51; 15 59; 16 13; 16 20; 16 22; 16 23; 16 36; 16 55; 17 3; 17 12; 17 20; 17 34; 17 41; 17 49; 17 56; 17 58; 18 8; 18 35; 18 41; 18 45; 19 2; 19 7; 19 24; 19 33; 19 49; 19 59; 20 6; 20 13; 20 22; 20 25; 20 27; 20 32; 20 38; 20 41; 20 46; 20 47; 20 57; 21 28; 21 33; 21 46; 21 57; 22 3; 22 5; 22 7; 22 11; 22 21; 22 24; 22 28; 22 33; 22 37; 22 40; 22 42; 22 58; 23 8; 23 21; 23 24; 23 44; 23 54; 23 57; 23 58; 24 14; 24 43; 24 45; 25 3; 25 7; 25 14; 25 22; 25 37; 25 41; 25 43; 25 52; 25 58; 26 5; 26 6; 26 48; 27 12; 27 39; 27 47; 28 30; 28 42; 28 46; 28 47; 29 26; 29 32; 29 53; 29 60; 30 19; 30 27; 30 38; 30 48; 30 52; 30 58; 31 10; 31 19; 31 22; 31 35; 31 44; 31 57; 32 17; 32 18; 32 43; 33 2; 33 4; 33 7; 33 10; 33 12; 33 36; 33 48; 33 59; 34 2; 34 10; 34 16; 34 38; 34 48; 34 57; 35 18; 35 22; 35 27; 35 37; 35 39; 35 57; 36 27; 36 34; 37 4; 37 49; 38 30; 38 43; 38 49; 39 3; 39 5; 40 3; 40 8; 40 15; 40 16; 40 17; 40 31; 40 41; 40 48; 40 55; 40 59; 41 27; 41 28; 41 30; 41 40; 42 18; 42 19; 42 27; 42 29; 42 33; 43 12; 43 33; 43 54; 43 56; 44 18; 44 28; 44 34; 44 40; 44 42; 44 52; 44 55; 45 17; 45 39; 45 40; 45 42; 45 52; 46 12; 46 13; 46 32; 46 35; 46 36; 46 47; 46 56; 46 57; 47 7; 47 13; 47 32; 47 45; 47 48; 47 53; 47 54; 47 57; 48 8; 48 14; 48 35; 48 53; 49 10; 49 21; 49 26; 50 5; 50 24; 50 28; 50 34; 50 36; 51 11; 51 17; 51 19; 51 25; 51 41; 51 52; 52 9; 52 14; 52 17; 52 20; 52 21; 52 34; 52 39; 52 51; 52 53; 52 59; 53 19; 53 28; 53 59; 54 11; 54 12; 54 29; 54 33; 54 41; 54 59; 55 16; 55 19; 55 40; 56 4; 56 5; 56 7; 56 10; 56 41; 56 46; 56 50; 57 16; 57 23; 57 31; 57 33; 57 34; 57 51; 58 3; 58 4; 58 17; 58 35; 58 47; 58 60; 59 9; 59 17; 59 20]
global d_x = [7.0, 7.0, 6.0, 5.0, 6.0, 3.0, 6.0, 2.0, 9.0, 1.0, 2.0, 8.0, 4.0, 8.0, 2.0, 5.0, 6.0, 8.0, 10.0, 8.0, 5.0, 8.0, 2.0, 10.0, 6.0, 8.0, 10.0, 7.0, 9.0, 8.0, 6.0, 3.0, 4.0, 4.0, 9.0, 10.0, 6.0, 9.0, 5.0, 10.0, 6.0, 6.0, 3.0, 1.0, 7.0, 9.0, 5.0, 3.0, 5.0, 8.0, 4.0, 4.0, 8.0, 1.0, 9.0, 8.0, 7.0, 7.0, 6.0, 10.0, 6.0, 6.0, 3.0, 4.0, 5.0, 5.0, 6.0, 5.0, 2.0, 6.0, 9.0, 8.0, 9.0, 4.0, 9.0, 2.0, 2.0, 5.0, 5.0, 8.0, 5.0, 8.0, 5.0, 5.0, 6.0, 7.0, 2.0, 2.0, 1.0, 10.0, 5.0, 1.0, 7.0, 10.0, 7.0, 2.0, 4.0, 5.0, 1.0, 9.0, 1.0, 8.0, 4.0, 2.0, 9.0, 10.0, 4.0, 1.0, 8.0, 10.0, 7.0, 7.0, 2.0, 4.0, 1.0, 8.0, 8.0, 10.0, 9.0, 10.0, 4.0, 10.0, 3.0, 1.0, 1.0, 4.0, 2.0, 3.0, 2.0, 3.0, 7.0, 4.0, 10.0, 7.0, 6.0, 7.0, 9.0, 8.0, 7.0, 8.0, 2.0, 1.0, 4.0, 8.0, 6.0, 5.0, 1.0, 7.0, 10.0, 10.0, 7.0, 4.0, 5.0, 8.0, 2.0, 3.0, 6.0, 6.0, 3.0, 3.0, 1.0, 10.0, 2.0, 6.0, 4.0, 3.0, 6.0, 5.0, 3.0, 4.0, 4.0, 4.0, 7.0, 6.0, 6.0, 7.0, 6.0, 5.0, 6.0, 9.0, 8.0, 3.0, 9.0, 9.0, 8.0, 1.0, 1.0, 4.0, 6.0, 3.0, 9.0, 9.0, 2.0, 1.0, 3.0, 6.0, 4.0, 6.0, 5.0, 1.0, 8.0, 3.0, 4.0, 8.0, 9.0, 2.0, 9.0, 9.0, 7.0, 1.0, 4.0, 2.0, 5.0, 5.0, 6.0, 4.0, 1.0, 8.0, 8.0, 5.0, 6.0, 3.0, 3.0, 7.0, 10.0, 6.0, 1.0, 6.0, 9.0, 1.0, 7.0, 7.0, 4.0, 10.0, 4.0, 3.0, 3.0, 2.0, 2.0, 8.0, 8.0, 10.0, 1.0, 9.0, 4.0, 4.0, 7.0, 6.0, 4.0, 9.0, 8.0, 10.0, 4.0, 6.0, 4.0, 3.0, 6.0, 9.0, 8.0, 9.0, 7.0, 9.0, 4.0, 3.0, 3.0, 3.0, 4.0, 2.0, 9.0, 2.0, 6.0, 10.0, 2.0, 1.0, 7.0, 6.0, 5.0, 10.0, 10.0, 3.0, 3.0, 6.0, 6.0, 10.0, 9.0, 10.0, 2.0, 9.0, 2.0, 3.0, 10.0, 5.0, 2.0, 3.0, 8.0, 3.0, 10.0, 9.0, 2.0, 3.0, 10.0, 3.0, 3.0, 2.0, 2.0, 4.0, 10.0, 2.0, 8.0, 7.0, 2.0, 10.0, 1.0, 9.0, 5.0, 1.0, 10.0, 5.0, 9.0, 10.0, 5.0, 5.0, 9.0, 6.0, 3.0, 4.0, 4.0, 10.0]
global b_x = 5
global d_y = [2.0, 9.0, 1.0, 3.0, 4.0, 8.0, 8.0, 5.0, 9.0, 3.0, 8.0, 4.0, 4.0, 6.0, 7.0, 6.0, 3.0, 6.0, 3.0, 4.0, 2.0, 2.0, 9.0, 4.0, 1.0, 3.0, 10.0, 6.0, 2.0, 5.0, 7.0, 10.0, 2.0, 2.0, 4.0, 2.0, 9.0, 6.0, 3.0, 4.0, 9.0, 10.0, 2.0, 4.0, 5.0, 8.0, 8.0, 4.0, 2.0, 4.0, 6.0, 7.0, 4.0, 7.0, 4.0, 8.0, 1.0, 2.0, 3.0, 9.0, 3.0, 5.0, 8.0, 5.0, 9.0, 7.0, 3.0, 6.0, 2.0, 7.0, 1.0, 1.0, 1.0, 6.0, 3.0, 1.0, 7.0, 10.0, 2.0, 10.0, 6.0, 6.0, 3.0, 5.0, 5.0, 9.0, 8.0, 2.0, 7.0, 2.0, 5.0, 8.0, 7.0, 10.0, 9.0, 10.0, 10.0, 8.0, 4.0, 6.0, 2.0, 8.0, 4.0, 3.0, 8.0, 4.0, 1.0, 9.0, 3.0, 8.0, 6.0, 1.0, 3.0, 7.0, 7.0, 10.0, 4.0, 8.0, 2.0, 6.0, 9.0, 6.0, 1.0, 3.0, 2.0, 2.0, 5.0, 3.0, 5.0, 1.0, 7.0, 2.0, 1.0, 2.0, 10.0, 3.0, 5.0, 5.0, 9.0, 8.0, 9.0, 4.0, 1.0, 9.0, 10.0, 2.0, 10.0, 5.0, 7.0, 5.0, 4.0, 10.0, 3.0, 9.0, 9.0, 1.0, 4.0, 6.0, 3.0, 7.0, 2.0, 7.0, 1.0, 7.0, 9.0, 5.0, 6.0, 1.0, 9.0, 4.0, 7.0, 9.0, 7.0, 9.0, 6.0, 7.0, 4.0, 5.0, 10.0, 3.0, 8.0, 5.0, 1.0, 10.0, 8.0, 10.0, 5.0, 7.0, 6.0, 6.0, 7.0, 1.0, 10.0, 2.0, 8.0, 2.0, 10.0, 7.0, 4.0, 6.0, 10.0, 7.0, 2.0, 8.0, 4.0, 3.0, 4.0, 7.0, 8.0, 1.0, 9.0, 2.0, 4.0, 10.0, 6.0, 3.0, 5.0, 1.0, 5.0, 10.0, 1.0, 10.0, 9.0, 5.0, 2.0, 8.0, 7.0, 3.0, 7.0, 2.0, 9.0, 10.0, 3.0, 1.0, 3.0, 5.0, 9.0, 10.0, 5.0, 5.0, 10.0, 6.0, 10.0, 8.0, 5.0, 3.0, 10.0, 3.0, 2.0, 4.0, 4.0, 5.0, 5.0, 8.0, 3.0, 7.0, 10.0, 9.0, 7.0, 5.0, 4.0, 1.0, 7.0, 7.0, 2.0, 5.0, 7.0, 1.0, 9.0, 1.0, 2.0, 4.0, 2.0, 8.0, 4.0, 2.0, 5.0, 2.0, 3.0, 4.0, 8.0, 5.0, 5.0, 8.0, 9.0, 8.0, 6.0, 8.0, 1.0, 5.0, 6.0, 6.0, 5.0, 9.0, 7.0, 9.0, 1.0, 9.0, 6.0, 7.0, 4.0, 10.0, 8.0, 3.0, 4.0, 4.0, 7.0, 10.0, 2.0, 4.0, 9.0, 2.0, 10.0, 4.0, 6.0, 7.0, 7.0, 6.0, 10.0, 3.0, 10.0, 4.0, 9.0, 10.0, 5.0, 9.0, 2.0, 7.0]
global b_y = 10
global p = [0.708, 0.693, 0.554, 0.158, 0.806, 0.497, 0.955, 0.553, 0.938, 0.495, 0.345, 0.492, 0.095, 0.791, 0.616, 0.673, 0.524, 0.996, 0.254, 0.159, 0.96, 0.511, 0.129, 0.741, 0.933, 0.111, 0.282, 0.354, 0.029, 0.898, 0.835, 0.239, 0.209, 0.146, 0.659, 0.597, 0.235, 0.692, 0.664, 0.652, 0.446, 0.577, 0.474, 0.616, 0.187, 0.341, 0.496, 0.492, 0.491, 0.617, 0.703, 0.713, 0.29, 0.448, 0.552, 0.32, 0.151, 0.398, 0.268, 0.74, 0.892, 0.553, 0.242, 0.928, 0.529, 0.203, 0.995, 0.468, 0.215, 0.39, 0.67, 0.648, 0.839, 0.294, 0.106, 0.177, 0.93, 0.267, 0.393, 0.681, 0.984, 0.866, 0.153, 0.357, 0.189, 0.35, 0.448, 0.351, 0.365, 0.306, 0.91, 0.03, 0.046, 0.093, 0.58, 0.772, 0.033, 0.426, 0.556, 0.319, 0.266, 0.459, 0.172, 0.49, 0.617, 0.378, 0.852, 0.105, 0.328, 0.118, 0.946, 0.561, 0.971, 0.112, 0.308, 0.292, 0.299, 0.911, 0.868, 0.639, 0.163, 0.857, 0.663, 0.094, 0.021, 0.061, 0.617, 0.467, 0.093, 0.539, 0.896, 0.899, 0.583, 0.113, 0.597, 0.037, 0.156, 0.352, 0.489, 0.374, 0.387, 0.128, 0.889, 0.154, 0.558, 0.39, 0.572, 0.666, 0.42, 0.671, 0.786, 0.587, 0.808, 0.905, 0.288, 0.973, 0.833, 0.315, 0.395, 0.973, 0.721, 0.123, 0.008, 0.18, 0.964, 0.725, 0.214, 0.932, 0.031, 0.764, 0.656, 0.252, 0.045, 0.158, 0.901, 0.256, 0.25, 0.448, 0.831, 0.649, 0.433, 0.554, 0.321, 0.402, 0.48, 0.226, 0.356, 0.012, 0.483, 0.692, 0.857, 0.706, 0.404, 0.569, 0.793, 0.512, 0.658, 0.04, 0.333, 0.625, 0.397, 0.648, 0.436, 0.451, 0.637, 0.035, 0.264, 0.184, 0.279, 0.921, 0.359, 0.397, 0.46, 0.917, 0.618, 0.392, 0.516, 0.806, 0.534, 0.397, 0.733, 0.714, 0.107, 0.797, 0.324, 0.451, 0.986, 0.703, 0.557, 0.457, 0.614, 0.14, 0.541, 0.966, 0.41, 0.54, 0.5, 0.456, 0.662, 0.837, 0.357, 0.537, 0.926, 0.546, 0.766, 0.53, 0.128, 0.19, 0.27, 0.56, 0.005, 0.19, 0.645, 0.278, 0.687, 0.512, 0.489, 0.991, 0.38, 0.168, 0.302, 0.083, 0.359, 0.742, 0.968, 0.481, 0.04, 0.961, 0.551, 0.672, 0.908, 0.224, 0.912, 0.225, 0.633, 0.176, 0.125, 0.903, 0.465, 0.877, 0.82, 0.265, 0.251, 0.616, 0.565, 0.764, 0.612, 0.531, 0.044, 0.885, 0.28, 0.777, 0.879, 0.261, 0.558, 0.593, 0.07, 0.13, 0.411, 0.57, 0.958, 0.917, 0.373, 0.063, 0.855, 0.845, 0.217, 0.888, 0.125, 0.27, 0.636, 0.349, 0.027, 0.02, 0.286, 0.421, 0.635, 0.752, 0.91, 0.257, 0.147, 0.561, 0.971, 0.827, 0.546, 0.322, 0.246, 0.28]
global q = [0.927, 0.791, 0.951, 0.329, 0.955, 0.512, 0.957, 0.812, 0.981, 0.783, 0.541, 0.621, 0.153, 0.946, 0.738, 0.917, 0.68, 0.999, 0.822, 0.842, 0.999, 0.678, 0.798, 0.951, 0.975, 0.157, 0.341, 0.875, 0.729, 0.961, 0.909, 0.734, 0.971, 0.992, 0.99, 0.718, 0.841, 0.765, 0.969, 0.874, 0.895, 0.679, 0.998, 0.754, 0.219, 0.682, 0.945, 0.624, 0.975, 0.861, 0.857, 0.986, 0.841, 0.665, 0.873, 0.959, 0.41, 0.755, 0.808, 0.838, 0.924, 0.778, 0.416, 0.966, 0.693, 0.274, 0.997, 0.978, 0.516, 0.505, 0.993, 0.722, 0.876, 0.357, 0.808, 0.824, 0.966, 0.371, 0.859, 0.794, 0.984, 0.952, 0.573, 0.596, 0.222, 0.462, 0.551, 0.696, 0.882, 0.481, 0.96, 0.523, 0.2, 0.763, 0.695, 0.81, 0.566, 0.877, 0.936, 0.806, 0.82, 0.459, 0.197, 0.685, 0.63, 0.857, 0.964, 0.544, 0.449, 0.786, 0.963, 0.75, 0.999, 0.604, 0.433, 0.304, 0.773, 0.992, 0.902, 0.662, 0.828, 0.926, 0.85, 0.676, 0.539, 0.92, 0.941, 0.618, 0.406, 0.882, 0.903, 0.972, 0.768, 0.324, 0.991, 0.176, 0.854, 0.469, 0.539, 0.464, 0.772, 0.576, 0.982, 0.472, 0.762, 0.578, 0.883, 0.795, 0.589, 0.818, 0.833, 0.795, 0.909, 0.933, 0.387, 0.981, 0.936, 0.384, 0.402, 0.99, 0.999, 0.951, 0.587, 0.257, 0.973, 0.924, 0.408, 0.966, 0.641, 0.995, 0.802, 0.355, 0.195, 0.257, 0.993, 0.692, 0.794, 0.772, 0.881, 0.822, 0.547, 0.649, 0.518, 0.707, 0.608, 0.57, 0.527, 0.093, 0.608, 0.878, 0.88, 0.765, 0.788, 0.957, 0.97, 0.731, 0.968, 0.349, 0.958, 0.936, 0.631, 0.696, 0.499, 0.818, 0.888, 0.463, 0.42, 0.363, 0.395, 0.953, 0.825, 0.805, 0.947, 0.985, 0.93, 0.398, 0.71, 0.828, 0.69, 0.692, 0.744, 0.839, 0.964, 0.997, 0.893, 0.889, 0.99, 0.815, 0.822, 0.839, 0.998, 0.77, 0.618, 0.996, 0.902, 0.56, 0.866, 0.635, 0.865, 0.942, 0.827, 0.887, 0.93, 0.638, 0.843, 0.91, 0.53, 0.267, 0.581, 0.883, 0.472, 0.684, 0.955, 0.643, 0.771, 0.98, 0.788, 0.996, 0.829, 0.903, 0.708, 0.454, 0.855, 0.946, 0.968, 0.504, 0.337, 0.972, 0.898, 0.829, 0.985, 0.319, 0.924, 0.791, 0.85, 0.616, 0.391, 0.99, 0.836, 0.934, 0.974, 0.724, 0.629, 0.669, 0.737, 0.935, 0.664, 0.786, 0.598, 0.997, 0.366, 0.884, 0.96, 0.694, 0.994, 0.599, 0.312, 0.421, 0.441, 0.913, 0.999, 0.979, 0.697, 0.395, 0.954, 0.894, 0.448, 0.922, 0.168, 0.825, 0.757, 0.612, 0.472, 0.085, 0.535, 0.467, 0.901, 0.77, 0.936, 0.491, 0.574, 0.57, 0.98, 0.922, 0.799, 0.858, 0.529, 0.796]
global origin = 1
global destination = 60