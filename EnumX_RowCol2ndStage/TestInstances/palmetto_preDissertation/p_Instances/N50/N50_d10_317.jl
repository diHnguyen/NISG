global arcs = [1 19; 1 26; 2 6; 2 10; 2 26; 2 32; 2 45; 3 12; 3 24; 3 35; 3 38; 4 2; 4 7; 5 10; 5 34; 5 45; 5 46; 6 17; 6 18; 6 21; 6 29; 6 31; 6 34; 6 45; 7 25; 7 27; 7 30; 7 35; 7 45; 7 46; 7 47; 8 11; 8 31; 8 40; 8 43; 8 45; 9 7; 9 28; 9 30; 9 41; 10 16; 10 32; 10 36; 10 43; 11 12; 11 19; 11 39; 11 40; 11 41; 12 10; 12 18; 12 20; 13 4; 13 11; 13 17; 13 26; 13 36; 13 38; 13 43; 13 49; 14 4; 14 28; 14 38; 14 39; 14 48; 15 3; 15 6; 15 8; 15 12; 15 19; 16 20; 16 23; 16 39; 16 40; 17 3; 17 20; 17 26; 17 31; 17 32; 18 10; 18 11; 18 30; 18 31; 18 49; 19 4; 19 14; 19 37; 20 10; 20 25; 20 27; 20 31; 20 35; 21 19; 21 26; 21 28; 22 2; 22 25; 22 32; 22 42; 23 20; 23 22; 23 44; 24 4; 24 20; 24 23; 24 29; 24 41; 24 43; 24 45; 24 48; 25 6; 25 7; 25 14; 25 30; 26 8; 26 28; 26 34; 26 50; 27 11; 27 22; 27 23; 27 38; 27 42; 27 44; 28 9; 28 23; 28 49; 29 3; 29 6; 29 23; 29 26; 29 42; 29 43; 30 7; 30 18; 30 41; 30 47; 31 7; 31 25; 31 26; 31 42; 31 49; 32 4; 32 9; 32 11; 32 13; 32 18; 32 20; 32 28; 32 29; 32 31; 32 38; 33 3; 33 4; 33 12; 33 21; 33 40; 34 8; 34 25; 34 32; 34 41; 35 12; 35 17; 35 22; 35 24; 35 36; 35 38; 35 41; 36 3; 36 33; 36 34; 36 35; 36 47; 37 4; 37 12; 37 24; 38 17; 38 29; 38 41; 38 46; 39 36; 39 40; 39 44; 40 3; 40 16; 40 24; 40 36; 40 44; 40 48; 41 27; 41 29; 41 33; 42 2; 42 7; 42 19; 42 26; 42 31; 43 5; 43 37; 43 45; 44 16; 45 17; 45 19; 45 40; 46 4; 46 8; 46 9; 46 13; 46 36; 46 39; 46 42; 46 47; 47 13; 47 15; 47 21; 47 23; 47 49; 47 50; 48 34; 49 2; 49 14; 49 21; 49 33; 49 38; 49 46; 49 47]
global d_x = [10.0, 10.0, 2.0, 1.0, 7.0, 3.0, 2.0, 6.0, 3.0, 7.0, 5.0, 10.0, 5.0, 5.0, 5.0, 3.0, 1.0, 6.0, 6.0, 1.0, 1.0, 6.0, 3.0, 6.0, 4.0, 8.0, 10.0, 7.0, 7.0, 5.0, 5.0, 2.0, 6.0, 2.0, 2.0, 8.0, 8.0, 2.0, 2.0, 8.0, 9.0, 8.0, 9.0, 9.0, 2.0, 10.0, 9.0, 10.0, 2.0, 6.0, 8.0, 5.0, 9.0, 2.0, 9.0, 2.0, 6.0, 3.0, 1.0, 1.0, 8.0, 10.0, 6.0, 8.0, 8.0, 10.0, 7.0, 7.0, 8.0, 10.0, 4.0, 2.0, 8.0, 7.0, 8.0, 9.0, 6.0, 2.0, 4.0, 2.0, 9.0, 7.0, 10.0, 6.0, 2.0, 1.0, 5.0, 10.0, 9.0, 2.0, 8.0, 8.0, 2.0, 7.0, 7.0, 6.0, 10.0, 6.0, 8.0, 10.0, 6.0, 2.0, 4.0, 9.0, 10.0, 10.0, 1.0, 3.0, 4.0, 10.0, 5.0, 9.0, 1.0, 9.0, 4.0, 10.0, 6.0, 10.0, 8.0, 3.0, 5.0, 5.0, 7.0, 7.0, 8.0, 6.0, 2.0, 8.0, 2.0, 7.0, 6.0, 9.0, 7.0, 7.0, 8.0, 4.0, 5.0, 6.0, 6.0, 9.0, 1.0, 1.0, 4.0, 4.0, 2.0, 1.0, 5.0, 4.0, 9.0, 6.0, 6.0, 5.0, 9.0, 10.0, 1.0, 9.0, 9.0, 10.0, 6.0, 1.0, 2.0, 10.0, 8.0, 6.0, 1.0, 5.0, 1.0, 7.0, 6.0, 10.0, 7.0, 7.0, 6.0, 4.0, 9.0, 8.0, 9.0, 8.0, 5.0, 7.0, 8.0, 6.0, 10.0, 2.0, 6.0, 7.0, 9.0, 5.0, 3.0, 9.0, 7.0, 5.0, 7.0, 5.0, 6.0, 6.0, 6.0, 4.0, 4.0, 7.0, 1.0, 1.0, 9.0, 6.0, 3.0, 1.0, 6.0, 9.0, 5.0, 3.0, 1.0, 10.0, 9.0, 4.0, 7.0, 8.0, 2.0, 7.0, 9.0, 5.0, 4.0, 7.0, 9.0, 6.0, 3.0, 6.0]
global b_x = 5
global d_y = [6.0, 5.0, 2.0, 7.0, 7.0, 10.0, 2.0, 9.0, 5.0, 6.0, 4.0, 3.0, 5.0, 3.0, 7.0, 2.0, 5.0, 9.0, 4.0, 6.0, 4.0, 8.0, 8.0, 4.0, 10.0, 6.0, 8.0, 1.0, 9.0, 4.0, 8.0, 6.0, 2.0, 6.0, 1.0, 8.0, 1.0, 3.0, 10.0, 1.0, 7.0, 1.0, 6.0, 5.0, 1.0, 4.0, 2.0, 3.0, 1.0, 9.0, 7.0, 5.0, 5.0, 1.0, 9.0, 4.0, 4.0, 5.0, 8.0, 2.0, 6.0, 1.0, 4.0, 7.0, 8.0, 9.0, 2.0, 9.0, 1.0, 5.0, 8.0, 3.0, 10.0, 6.0, 4.0, 5.0, 3.0, 3.0, 7.0, 10.0, 1.0, 4.0, 6.0, 6.0, 2.0, 4.0, 8.0, 7.0, 10.0, 10.0, 4.0, 2.0, 2.0, 10.0, 4.0, 10.0, 3.0, 1.0, 8.0, 6.0, 9.0, 9.0, 2.0, 10.0, 10.0, 5.0, 8.0, 1.0, 9.0, 2.0, 10.0, 6.0, 10.0, 5.0, 7.0, 1.0, 1.0, 9.0, 2.0, 10.0, 1.0, 2.0, 8.0, 9.0, 5.0, 1.0, 1.0, 1.0, 2.0, 8.0, 1.0, 10.0, 4.0, 1.0, 6.0, 1.0, 2.0, 5.0, 7.0, 5.0, 5.0, 9.0, 6.0, 6.0, 8.0, 8.0, 7.0, 6.0, 5.0, 6.0, 1.0, 10.0, 3.0, 10.0, 10.0, 8.0, 3.0, 8.0, 4.0, 7.0, 4.0, 8.0, 10.0, 9.0, 7.0, 2.0, 8.0, 2.0, 7.0, 2.0, 7.0, 3.0, 10.0, 4.0, 3.0, 1.0, 2.0, 7.0, 3.0, 7.0, 10.0, 6.0, 9.0, 8.0, 10.0, 6.0, 8.0, 6.0, 6.0, 9.0, 9.0, 8.0, 2.0, 5.0, 7.0, 2.0, 1.0, 10.0, 2.0, 10.0, 7.0, 9.0, 5.0, 1.0, 1.0, 10.0, 6.0, 9.0, 2.0, 3.0, 5.0, 4.0, 4.0, 5.0, 10.0, 1.0, 7.0, 5.0, 5.0, 2.0, 10.0, 6.0, 2.0, 8.0, 10.0, 4.0]
global b_y = 10
global p = [0.444, 0.9, 0.227, 0.617, 0.6, 0.502, 0.081, 0.374, 0.152, 0.209, 0.436, 0.784, 0.94, 0.347, 0.086, 0.106, 0.728, 0.276, 0.219, 0.696, 0.639, 0.28, 0.174, 0.555, 0.051, 0.858, 0.569, 0.588, 0.939, 0.218, 0.858, 0.821, 0.047, 0.652, 0.006, 0.153, 0.377, 0.47, 0.118, 0.84, 0.055, 0.494, 0.763, 0.959, 0.726, 0.483, 0.083, 0.319, 0.948, 0.888, 0.117, 0.028, 0.982, 0.86, 0.521, 0.168, 0.968, 0.541, 0.066, 0.736, 0.302, 0.892, 0.047, 0.274, 0.355, 0.487, 0.574, 0.721, 0.106, 0.641, 0.071, 0.193, 0.198, 0.991, 0.102, 0.289, 0.155, 0.691, 0.303, 0.269, 0.747, 0.309, 0.899, 0.672, 0.677, 0.271, 0.712, 0.121, 0.09, 0.761, 0.953, 0.911, 0.954, 0.094, 0.831, 0.859, 0.436, 0.214, 0.769, 0.704, 0.373, 0.591, 0.555, 0.283, 0.785, 0.627, 0.462, 0.71, 0.239, 0.712, 0.397, 0.619, 0.395, 0.778, 0.272, 0.362, 0.273, 0.89, 0.721, 0.1, 0.144, 0.512, 0.421, 0.089, 0.636, 0.782, 0.31, 0.333, 0.063, 0.592, 0.495, 0.432, 0.476, 0.956, 0.979, 0.737, 0.528, 0.782, 0.621, 0.597, 0.888, 0.465, 0.488, 0.231, 0.429, 0.938, 0.471, 0.916, 0.974, 0.671, 0.167, 0.171, 0.509, 0.27, 0.511, 0.417, 0.947, 0.204, 0.79, 0.976, 0.668, 0.509, 0.868, 0.433, 0.148, 0.531, 0.15, 0.738, 0.147, 0.251, 0.716, 0.637, 0.004, 0.6, 0.741, 0.214, 0.628, 0.157, 0.538, 0.61, 0.594, 0.904, 0.926, 0.7, 0.532, 0.398, 0.785, 0.176, 0.16, 0.679, 0.56, 0.382, 0.9, 0.059, 0.146, 0.604, 0.429, 0.32, 0.949, 0.802, 0.782, 0.542, 0.208, 0.349, 0.091, 0.605, 0.564, 0.465, 0.063, 0.718, 0.607, 0.299, 0.678, 0.083, 0.479, 0.048, 0.105, 0.618, 0.589, 0.489, 0.637, 0.13, 0.395, 0.462, 0.706, 0.954]
global q = [0.535, 0.933, 0.415, 0.804, 0.625, 0.634, 0.62, 0.462, 0.886, 0.35, 0.974, 0.844, 0.95, 0.85, 0.824, 0.523, 0.921, 0.897, 0.428, 0.708, 0.82, 0.98, 0.455, 0.802, 0.707, 0.994, 0.824, 0.731, 0.946, 0.743, 0.968, 0.864, 0.049, 0.845, 0.726, 0.685, 0.387, 0.841, 0.539, 0.902, 0.332, 0.524, 0.995, 0.965, 0.933, 0.514, 0.724, 0.92, 0.998, 0.9, 0.771, 0.351, 0.994, 0.941, 0.663, 0.369, 0.996, 0.69, 0.747, 0.771, 0.599, 0.917, 0.105, 0.277, 0.613, 0.889, 0.885, 0.997, 0.352, 0.937, 0.243, 0.967, 0.683, 0.996, 0.622, 0.788, 0.622, 0.972, 0.767, 0.533, 0.776, 0.387, 0.981, 0.681, 0.943, 0.385, 0.795, 0.324, 0.511, 0.844, 0.958, 0.912, 0.988, 0.963, 0.84, 0.978, 0.862, 0.325, 0.916, 0.843, 0.454, 0.755, 0.951, 0.633, 0.866, 0.63, 0.921, 0.92, 0.937, 0.783, 0.423, 0.838, 0.711, 0.945, 0.888, 0.539, 0.642, 0.893, 0.833, 0.16, 0.172, 0.76, 0.817, 0.191, 0.852, 0.964, 0.382, 0.769, 0.961, 0.869, 0.735, 0.853, 0.537, 0.974, 0.988, 0.916, 0.766, 0.849, 0.824, 0.714, 0.907, 0.781, 0.69, 0.231, 0.875, 0.971, 0.961, 0.94, 0.986, 0.922, 0.961, 0.949, 0.881, 0.287, 0.517, 0.468, 0.967, 0.692, 0.867, 0.989, 0.872, 0.861, 0.988, 0.715, 0.753, 0.769, 0.202, 0.989, 0.591, 0.505, 0.953, 0.998, 0.978, 0.804, 0.816, 0.9, 0.775, 0.321, 0.579, 0.803, 0.606, 0.969, 0.96, 0.798, 0.929, 0.926, 0.991, 0.549, 0.615, 0.824, 0.78, 0.695, 0.95, 0.826, 0.362, 0.963, 0.888, 0.337, 0.967, 0.896, 0.852, 0.849, 0.565, 0.966, 0.558, 0.664, 0.99, 0.548, 0.674, 0.817, 0.647, 0.321, 0.877, 0.929, 0.574, 0.725, 0.156, 0.646, 0.619, 0.712, 0.915, 0.275, 0.904, 0.797, 0.768, 0.97]
global origin = 1
global destination = 50