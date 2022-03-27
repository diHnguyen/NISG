global arcs = [1 6; 1 10; 1 17; 1 21; 1 23; 1 27; 1 32; 1 43; 2 5; 2 23; 2 33; 2 39; 2 46; 3 10; 3 12; 3 18; 3 24; 3 29; 3 34; 3 42; 3 48; 4 7; 4 10; 4 13; 4 20; 4 38; 5 2; 5 11; 5 18; 5 42; 6 4; 6 5; 6 8; 6 22; 6 32; 6 35; 6 36; 6 38; 7 16; 7 17; 7 23; 7 29; 7 33; 7 48; 8 11; 8 35; 8 39; 8 40; 8 45; 8 48; 9 8; 9 21; 9 28; 9 31; 9 39; 9 50; 10 6; 10 32; 11 13; 11 15; 11 26; 11 39; 11 46; 12 6; 12 8; 12 13; 12 19; 12 36; 12 38; 13 22; 13 26; 13 28; 13 41; 13 42; 14 4; 14 42; 14 44; 14 45; 14 46; 14 48; 14 50; 15 7; 15 24; 15 30; 15 33; 15 49; 16 9; 16 14; 16 32; 16 45; 17 24; 17 39; 18 2; 18 4; 18 12; 18 16; 18 21; 18 30; 18 35; 19 21; 19 23; 19 33; 20 5; 20 13; 20 28; 20 30; 20 40; 21 7; 21 33; 21 46; 22 16; 23 25; 23 28; 23 29; 23 36; 24 4; 24 7; 24 15; 24 41; 25 10; 25 11; 25 40; 25 47; 25 50; 26 3; 26 6; 26 16; 26 27; 26 32; 26 40; 26 46; 27 4; 27 12; 27 26; 27 42; 28 2; 28 13; 28 14; 28 45; 28 50; 29 12; 29 14; 29 48; 30 3; 30 5; 30 7; 30 21; 30 39; 30 46; 30 48; 31 4; 31 8; 31 23; 32 9; 32 11; 32 16; 32 29; 32 33; 32 35; 32 46; 33 43; 33 44; 34 2; 34 48; 35 6; 35 11; 35 12; 35 33; 35 42; 35 43; 36 23; 36 25; 36 30; 36 39; 37 8; 37 23; 37 42; 38 5; 38 13; 38 18; 38 19; 38 21; 38 31; 38 39; 38 43; 38 44; 39 4; 39 11; 39 16; 39 17; 39 34; 39 43; 40 20; 40 22; 40 25; 40 32; 40 37; 40 43; 40 46; 41 5; 41 29; 41 36; 41 49; 42 22; 42 24; 42 25; 42 32; 42 46; 42 49; 43 4; 43 12; 43 13; 43 15; 44 3; 44 5; 44 25; 44 41; 44 48; 45 2; 45 21; 45 26; 45 33; 45 34; 46 4; 46 6; 46 8; 46 14; 46 16; 46 18; 46 28; 46 34; 46 37; 47 7; 47 26; 48 6; 48 16; 48 22; 48 23; 48 29; 48 36; 48 40; 48 46; 49 7; 49 18; 49 26; 49 27; 49 31; 49 33]
global d_x = [5.0, 5.0, 9.0, 6.0, 9.0, 9.0, 7.0, 6.0, 6.0, 6.0, 2.0, 9.0, 6.0, 2.0, 4.0, 5.0, 10.0, 4.0, 3.0, 6.0, 1.0, 8.0, 3.0, 9.0, 4.0, 6.0, 10.0, 10.0, 2.0, 2.0, 8.0, 3.0, 2.0, 9.0, 9.0, 6.0, 3.0, 7.0, 3.0, 10.0, 4.0, 4.0, 1.0, 1.0, 3.0, 9.0, 10.0, 5.0, 8.0, 8.0, 3.0, 8.0, 2.0, 6.0, 9.0, 2.0, 7.0, 10.0, 2.0, 10.0, 9.0, 10.0, 3.0, 4.0, 6.0, 2.0, 7.0, 7.0, 9.0, 7.0, 3.0, 9.0, 10.0, 9.0, 5.0, 4.0, 10.0, 3.0, 5.0, 7.0, 2.0, 10.0, 1.0, 8.0, 7.0, 8.0, 3.0, 6.0, 7.0, 6.0, 8.0, 4.0, 3.0, 8.0, 10.0, 3.0, 5.0, 5.0, 9.0, 1.0, 2.0, 9.0, 4.0, 7.0, 2.0, 6.0, 9.0, 1.0, 9.0, 1.0, 8.0, 2.0, 10.0, 10.0, 3.0, 7.0, 2.0, 1.0, 4.0, 6.0, 1.0, 10.0, 1.0, 10.0, 1.0, 3.0, 7.0, 3.0, 10.0, 7.0, 4.0, 5.0, 1.0, 9.0, 7.0, 7.0, 6.0, 10.0, 9.0, 9.0, 6.0, 6.0, 3.0, 2.0, 6.0, 2.0, 4.0, 3.0, 6.0, 2.0, 10.0, 2.0, 8.0, 8.0, 9.0, 8.0, 9.0, 7.0, 7.0, 10.0, 10.0, 10.0, 7.0, 10.0, 8.0, 1.0, 8.0, 7.0, 5.0, 1.0, 10.0, 6.0, 4.0, 3.0, 5.0, 3.0, 7.0, 8.0, 5.0, 5.0, 2.0, 10.0, 2.0, 3.0, 3.0, 4.0, 4.0, 9.0, 9.0, 3.0, 10.0, 1.0, 4.0, 5.0, 6.0, 3.0, 7.0, 8.0, 4.0, 2.0, 10.0, 3.0, 5.0, 4.0, 1.0, 1.0, 2.0, 5.0, 10.0, 8.0, 10.0, 10.0, 9.0, 6.0, 5.0, 7.0, 2.0, 10.0, 5.0, 10.0, 5.0, 3.0, 10.0, 7.0, 4.0, 3.0, 8.0, 3.0, 1.0, 2.0, 2.0, 1.0, 7.0, 5.0, 5.0, 4.0, 5.0, 10.0, 2.0, 5.0, 10.0, 3.0, 8.0, 5.0, 10.0, 2.0, 3.0, 4.0]
global b_x = 5
global d_y = [7.0, 3.0, 2.0, 4.0, 8.0, 10.0, 6.0, 2.0, 3.0, 1.0, 9.0, 1.0, 8.0, 9.0, 10.0, 7.0, 4.0, 4.0, 2.0, 3.0, 9.0, 5.0, 8.0, 5.0, 10.0, 2.0, 6.0, 3.0, 9.0, 10.0, 8.0, 6.0, 1.0, 2.0, 3.0, 7.0, 8.0, 1.0, 1.0, 5.0, 4.0, 7.0, 6.0, 7.0, 1.0, 5.0, 8.0, 2.0, 1.0, 5.0, 3.0, 3.0, 6.0, 4.0, 5.0, 8.0, 5.0, 4.0, 2.0, 3.0, 2.0, 7.0, 10.0, 4.0, 1.0, 4.0, 6.0, 6.0, 2.0, 9.0, 8.0, 3.0, 8.0, 6.0, 3.0, 3.0, 9.0, 6.0, 9.0, 1.0, 7.0, 9.0, 7.0, 3.0, 5.0, 10.0, 8.0, 7.0, 5.0, 2.0, 5.0, 8.0, 9.0, 8.0, 3.0, 7.0, 10.0, 4.0, 4.0, 10.0, 7.0, 4.0, 7.0, 6.0, 5.0, 3.0, 8.0, 6.0, 3.0, 8.0, 4.0, 9.0, 7.0, 3.0, 5.0, 9.0, 4.0, 10.0, 9.0, 7.0, 4.0, 8.0, 1.0, 8.0, 9.0, 5.0, 1.0, 10.0, 5.0, 9.0, 7.0, 10.0, 6.0, 7.0, 3.0, 5.0, 6.0, 1.0, 6.0, 1.0, 8.0, 2.0, 3.0, 2.0, 7.0, 2.0, 7.0, 10.0, 3.0, 7.0, 5.0, 7.0, 2.0, 9.0, 10.0, 8.0, 4.0, 5.0, 1.0, 2.0, 10.0, 7.0, 1.0, 5.0, 5.0, 6.0, 8.0, 4.0, 2.0, 7.0, 8.0, 5.0, 8.0, 9.0, 9.0, 6.0, 4.0, 5.0, 8.0, 1.0, 10.0, 3.0, 7.0, 9.0, 3.0, 5.0, 1.0, 6.0, 4.0, 6.0, 9.0, 8.0, 8.0, 4.0, 10.0, 4.0, 6.0, 3.0, 8.0, 8.0, 4.0, 7.0, 6.0, 4.0, 8.0, 4.0, 8.0, 8.0, 3.0, 7.0, 6.0, 5.0, 2.0, 3.0, 10.0, 5.0, 7.0, 8.0, 9.0, 9.0, 2.0, 5.0, 7.0, 4.0, 10.0, 2.0, 6.0, 3.0, 6.0, 9.0, 7.0, 7.0, 5.0, 8.0, 4.0, 5.0, 8.0, 3.0, 6.0, 7.0, 3.0, 6.0, 1.0, 2.0, 1.0, 8.0, 7.0, 5.0]
global b_y = 10
global p = [0.553, 0.151, 0.117, 0.84, 0.514, 0.943, 0.36, 0.707, 0.263, 0.112, 0.696, 0.04, 0.561, 0.955, 0.892, 0.919, 0.697, 0.956, 0.694, 0.303, 0.844, 0.381, 0.454, 0.394, 0.902, 0.577, 0.31, 0.587, 0.858, 0.742, 0.376, 0.864, 0.963, 0.914, 0.295, 0.237, 0.82, 0.233, 0.362, 0.44, 0.147, 0.092, 0.925, 0.21, 0.127, 0.338, 0.994, 0.271, 0.553, 0.823, 0.699, 0.285, 0.796, 0.415, 0.859, 0.99, 0.093, 0.443, 0.097, 0.277, 0.977, 0.718, 0.239, 0.292, 0.172, 0.814, 0.338, 0.019, 0.633, 0.447, 0.727, 0.849, 0.408, 0.884, 0.488, 0.372, 0.66, 0.255, 0.089, 0.358, 0.392, 0.688, 0.184, 0.737, 0.047, 0.873, 0.61, 0.988, 0.177, 0.347, 0.588, 0.407, 0.813, 0.282, 0.869, 0.445, 0.887, 0.471, 0.092, 0.809, 0.959, 0.348, 0.559, 0.296, 0.737, 0.37, 0.406, 0.345, 0.242, 0.016, 0.9, 0.971, 0.052, 0.104, 0.638, 0.353, 0.586, 0.108, 0.984, 0.152, 0.602, 0.361, 0.684, 0.44, 0.144, 0.433, 0.813, 0.284, 0.173, 0.494, 0.224, 0.063, 0.316, 0.116, 0.883, 0.969, 0.407, 0.386, 0.435, 0.092, 0.634, 0.46, 0.605, 0.987, 0.24, 0.744, 0.53, 0.311, 0.632, 0.575, 0.453, 0.061, 0.055, 0.03, 0.485, 0.241, 0.91, 0.403, 0.542, 0.141, 0.428, 0.739, 0.662, 0.445, 0.367, 0.156, 0.363, 0.119, 0.607, 0.313, 0.77, 0.86, 0.923, 0.746, 0.963, 0.415, 0.554, 0.257, 0.889, 0.11, 0.468, 0.252, 0.031, 0.274, 0.079, 0.042, 0.354, 0.066, 0.029, 0.793, 0.124, 0.095, 0.744, 0.173, 0.715, 0.242, 0.233, 0.34, 0.81, 0.913, 0.741, 0.214, 0.997, 0.601, 0.975, 0.268, 0.068, 0.896, 0.972, 0.167, 0.452, 0.028, 0.292, 0.202, 0.778, 0.812, 0.432, 0.531, 0.485, 0.05, 0.965, 0.717, 0.389, 0.069, 0.746, 0.488, 0.255, 0.273, 0.438, 0.624, 0.269, 0.67, 0.784, 0.093, 0.966, 0.285, 0.636, 0.596, 0.933, 0.329, 0.067, 0.01, 0.447, 0.98, 0.091, 0.572, 0.675, 0.316]
global q = [0.583, 0.41, 0.643, 0.856, 0.587, 0.952, 0.538, 0.74, 0.358, 0.206, 0.881, 0.48, 0.955, 0.967, 0.961, 0.942, 0.883, 0.989, 0.783, 0.601, 0.886, 0.983, 0.557, 0.653, 0.954, 0.969, 0.968, 0.822, 0.865, 0.767, 0.622, 0.918, 0.993, 0.963, 0.604, 0.853, 0.94, 0.396, 0.476, 0.83, 0.264, 0.778, 0.973, 0.34, 0.204, 0.402, 0.998, 0.917, 0.644, 0.925, 0.814, 0.51, 0.904, 0.971, 0.908, 0.997, 0.166, 0.839, 0.664, 0.534, 0.996, 0.726, 0.56, 0.882, 0.258, 0.915, 0.938, 0.884, 0.868, 0.505, 0.893, 0.865, 0.624, 0.973, 0.922, 0.803, 0.967, 0.382, 0.118, 0.54, 0.494, 0.972, 0.532, 0.872, 0.684, 0.96, 0.955, 0.99, 0.732, 0.378, 0.794, 0.767, 0.851, 0.662, 0.877, 0.649, 0.958, 0.882, 0.105, 0.931, 0.977, 0.418, 0.862, 0.439, 0.789, 0.699, 0.408, 0.463, 0.257, 0.262, 0.955, 0.981, 0.055, 0.748, 0.675, 0.441, 0.645, 0.575, 0.995, 0.821, 0.66, 0.583, 0.911, 0.868, 0.678, 0.953, 0.873, 0.965, 0.329, 0.681, 0.471, 0.557, 0.724, 0.224, 0.998, 0.996, 0.854, 0.711, 0.941, 0.578, 0.735, 0.815, 0.757, 0.993, 0.486, 0.969, 0.63, 0.467, 0.742, 0.708, 0.486, 0.805, 0.532, 0.282, 0.736, 0.649, 0.951, 0.712, 0.716, 0.865, 0.982, 0.908, 0.75, 0.495, 0.806, 0.614, 0.582, 0.814, 0.71, 0.862, 0.854, 0.888, 0.964, 0.811, 0.997, 0.422, 0.562, 0.616, 0.986, 0.56, 0.655, 0.564, 0.221, 0.445, 0.917, 0.985, 0.52, 0.663, 0.778, 0.919, 0.478, 0.82, 0.804, 0.22, 0.736, 0.516, 0.325, 0.39, 0.811, 0.987, 0.741, 0.953, 0.997, 0.999, 0.985, 0.4, 0.113, 0.925, 0.98, 0.9, 0.836, 0.223, 0.322, 0.892, 0.921, 0.972, 0.797, 0.876, 0.524, 0.971, 0.993, 0.798, 0.931, 0.073, 0.865, 0.614, 0.818, 0.768, 0.803, 0.953, 0.991, 0.893, 0.916, 0.904, 0.991, 0.296, 0.8, 0.767, 0.996, 0.889, 0.964, 0.927, 0.991, 0.997, 0.262, 0.717, 0.73, 0.463]
global origin = 1
global destination = 50