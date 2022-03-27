global arcs = [1 3; 1 28; 1 33; 1 39; 2 6; 2 11; 2 20; 2 26; 2 30; 2 39; 3 10; 3 18; 3 22; 3 29; 3 50; 4 12; 4 20; 4 24; 4 38; 4 44; 5 4; 5 10; 5 27; 5 35; 6 2; 6 28; 6 40; 7 12; 7 13; 7 14; 7 20; 7 33; 7 48; 8 12; 8 24; 8 35; 8 37; 8 41; 8 45; 8 47; 9 28; 9 35; 9 40; 9 46; 10 3; 10 15; 10 21; 11 8; 11 10; 11 13; 11 26; 11 30; 12 8; 12 18; 12 26; 12 42; 13 4; 13 9; 13 37; 13 43; 13 49; 14 3; 14 16; 14 22; 14 34; 14 35; 14 48; 15 18; 15 22; 15 25; 15 26; 15 34; 15 42; 15 43; 15 46; 15 48; 16 45; 17 5; 17 6; 17 32; 17 37; 17 41; 17 43; 17 46; 18 41; 18 48; 19 10; 19 18; 20 15; 20 39; 20 47; 21 38; 21 43; 22 5; 22 13; 22 18; 23 6; 23 15; 23 22; 23 26; 23 27; 23 34; 23 46; 23 49; 24 5; 24 14; 24 27; 25 13; 25 26; 25 40; 25 42; 26 10; 26 11; 26 13; 26 15; 26 25; 26 48; 27 14; 27 16; 27 21; 27 34; 27 45; 28 12; 28 14; 28 31; 28 45; 29 5; 29 7; 29 10; 29 19; 29 34; 29 37; 29 39; 29 50; 30 3; 30 6; 30 35; 30 46; 31 10; 31 11; 31 21; 31 24; 31 32; 31 44; 32 7; 32 8; 32 25; 32 35; 33 4; 33 9; 33 42; 34 48; 35 11; 35 22; 35 25; 35 28; 35 49; 36 8; 36 11; 36 13; 36 20; 36 29; 36 30; 36 32; 36 49; 37 17; 37 26; 37 46; 38 3; 38 7; 38 20; 38 23; 38 30; 38 43; 38 48; 39 8; 39 25; 39 28; 39 31; 40 8; 40 19; 40 33; 40 48; 41 19; 41 21; 41 29; 42 22; 42 25; 42 26; 42 27; 42 32; 42 45; 43 6; 43 17; 43 28; 43 31; 43 35; 43 42; 43 50; 44 13; 44 19; 44 27; 44 36; 45 20; 45 33; 46 5; 46 17; 46 40; 47 4; 47 9; 47 20; 47 29; 47 32; 48 10; 48 12; 48 15; 48 22; 48 35; 49 3; 49 4; 49 8; 49 9; 49 12; 49 28; 49 29]
global d_x = [8.0, 3.0, 6.0, 2.0, 7.0, 1.0, 10.0, 1.0, 8.0, 2.0, 3.0, 1.0, 8.0, 2.0, 7.0, 5.0, 4.0, 4.0, 6.0, 1.0, 7.0, 7.0, 1.0, 6.0, 10.0, 2.0, 5.0, 9.0, 6.0, 1.0, 6.0, 6.0, 7.0, 1.0, 7.0, 4.0, 1.0, 7.0, 2.0, 8.0, 7.0, 10.0, 10.0, 7.0, 5.0, 6.0, 6.0, 4.0, 7.0, 1.0, 6.0, 3.0, 8.0, 1.0, 9.0, 6.0, 3.0, 3.0, 1.0, 3.0, 10.0, 7.0, 9.0, 4.0, 4.0, 6.0, 7.0, 8.0, 5.0, 3.0, 8.0, 7.0, 2.0, 2.0, 9.0, 9.0, 6.0, 6.0, 4.0, 5.0, 4.0, 1.0, 8.0, 7.0, 3.0, 3.0, 7.0, 9.0, 6.0, 1.0, 4.0, 9.0, 10.0, 7.0, 2.0, 10.0, 9.0, 1.0, 4.0, 4.0, 6.0, 8.0, 8.0, 10.0, 3.0, 3.0, 3.0, 4.0, 10.0, 2.0, 6.0, 5.0, 10.0, 10.0, 10.0, 3.0, 6.0, 5.0, 8.0, 10.0, 8.0, 4.0, 9.0, 1.0, 7.0, 2.0, 2.0, 4.0, 4.0, 8.0, 3.0, 6.0, 10.0, 3.0, 5.0, 9.0, 8.0, 7.0, 3.0, 10.0, 2.0, 8.0, 1.0, 1.0, 9.0, 5.0, 10.0, 7.0, 10.0, 3.0, 9.0, 7.0, 3.0, 5.0, 3.0, 1.0, 5.0, 6.0, 6.0, 3.0, 10.0, 10.0, 3.0, 7.0, 3.0, 7.0, 8.0, 5.0, 4.0, 8.0, 4.0, 9.0, 1.0, 9.0, 4.0, 4.0, 6.0, 6.0, 9.0, 4.0, 2.0, 7.0, 1.0, 10.0, 6.0, 9.0, 3.0, 10.0, 2.0, 7.0, 4.0, 2.0, 8.0, 9.0, 9.0, 10.0, 5.0, 3.0, 5.0, 6.0, 8.0, 6.0, 10.0, 5.0, 2.0, 2.0, 9.0, 4.0, 5.0, 6.0, 1.0, 8.0, 3.0, 7.0, 3.0, 3.0, 10.0, 6.0, 3.0, 8.0, 2.0, 9.0, 10.0, 6.0, 4.0]
global b_x = 5
global d_y = [8.0, 5.0, 5.0, 10.0, 8.0, 5.0, 8.0, 10.0, 1.0, 5.0, 7.0, 2.0, 5.0, 4.0, 2.0, 2.0, 6.0, 2.0, 6.0, 4.0, 1.0, 4.0, 4.0, 9.0, 1.0, 3.0, 9.0, 1.0, 2.0, 1.0, 2.0, 6.0, 8.0, 5.0, 6.0, 7.0, 4.0, 5.0, 5.0, 7.0, 3.0, 1.0, 10.0, 2.0, 8.0, 8.0, 2.0, 8.0, 10.0, 2.0, 4.0, 6.0, 4.0, 9.0, 3.0, 8.0, 3.0, 3.0, 4.0, 1.0, 7.0, 3.0, 8.0, 10.0, 3.0, 4.0, 5.0, 5.0, 7.0, 3.0, 5.0, 7.0, 10.0, 1.0, 8.0, 6.0, 10.0, 9.0, 2.0, 9.0, 6.0, 8.0, 1.0, 5.0, 4.0, 9.0, 6.0, 1.0, 1.0, 9.0, 5.0, 9.0, 3.0, 1.0, 4.0, 6.0, 5.0, 7.0, 4.0, 3.0, 8.0, 6.0, 3.0, 7.0, 5.0, 6.0, 6.0, 3.0, 4.0, 2.0, 9.0, 3.0, 10.0, 4.0, 7.0, 5.0, 7.0, 2.0, 9.0, 4.0, 2.0, 7.0, 10.0, 2.0, 5.0, 6.0, 10.0, 10.0, 2.0, 1.0, 2.0, 9.0, 1.0, 10.0, 3.0, 9.0, 3.0, 6.0, 10.0, 9.0, 3.0, 7.0, 3.0, 9.0, 4.0, 4.0, 3.0, 8.0, 3.0, 1.0, 6.0, 7.0, 3.0, 1.0, 10.0, 1.0, 7.0, 1.0, 6.0, 4.0, 8.0, 3.0, 8.0, 4.0, 6.0, 7.0, 2.0, 7.0, 5.0, 5.0, 10.0, 8.0, 7.0, 10.0, 5.0, 2.0, 6.0, 9.0, 2.0, 3.0, 1.0, 10.0, 6.0, 6.0, 2.0, 10.0, 9.0, 3.0, 1.0, 1.0, 5.0, 5.0, 1.0, 10.0, 6.0, 3.0, 5.0, 4.0, 4.0, 2.0, 9.0, 3.0, 1.0, 9.0, 6.0, 9.0, 10.0, 3.0, 4.0, 3.0, 3.0, 4.0, 8.0, 3.0, 7.0, 2.0, 4.0, 4.0, 9.0, 3.0, 10.0, 7.0, 6.0, 6.0, 3.0]
global b_y = 10
global p = [0.019, 0.707, 0.621, 0.163, 0.545, 0.127, 0.703, 0.807, 0.072, 0.168, 0.933, 0.711, 0.294, 0.37, 0.086, 0.132, 0.645, 0.211, 0.134, 0.851, 0.297, 0.851, 0.076, 0.95, 0.044, 0.224, 0.494, 0.159, 0.468, 0.601, 0.793, 0.509, 0.347, 0.531, 0.27, 0.738, 0.142, 0.171, 0.887, 0.343, 0.961, 0.838, 0.047, 0.622, 0.08, 0.697, 0.804, 0.048, 0.586, 0.715, 0.912, 0.43, 0.992, 0.352, 0.057, 0.155, 0.071, 0.408, 0.664, 0.087, 0.541, 0.617, 0.536, 0.069, 0.733, 0.518, 0.399, 0.351, 0.976, 0.628, 0.79, 0.238, 0.689, 0.489, 0.49, 0.012, 0.4, 0.505, 0.2, 0.28, 0.463, 0.033, 0.598, 0.888, 0.189, 0.767, 0.704, 0.025, 0.871, 0.182, 0.646, 0.55, 0.188, 0.317, 0.95, 0.847, 0.185, 0.602, 0.299, 0.808, 0.117, 0.242, 0.625, 0.195, 0.537, 0.66, 0.844, 0.31, 0.185, 0.598, 0.897, 0.852, 0.815, 0.413, 0.323, 0.251, 0.952, 0.237, 0.793, 0.183, 0.748, 0.158, 0.231, 0.962, 0.04, 0.33, 0.026, 0.979, 0.281, 0.137, 0.182, 0.918, 0.075, 0.58, 0.033, 0.39, 0.799, 0.299, 0.441, 0.817, 0.538, 0.651, 0.698, 0.837, 0.638, 0.649, 0.353, 0.377, 0.712, 0.297, 0.697, 0.552, 0.35, 0.202, 0.211, 0.167, 0.643, 0.879, 0.738, 0.456, 0.198, 0.558, 0.704, 0.293, 0.818, 0.186, 0.822, 0.601, 0.313, 0.179, 0.597, 0.065, 0.329, 0.033, 0.167, 0.648, 0.001, 0.388, 0.736, 0.802, 0.939, 0.598, 0.53, 0.416, 0.126, 0.196, 0.007, 0.706, 0.092, 0.469, 0.755, 0.019, 0.33, 0.343, 0.273, 0.132, 0.047, 0.197, 0.296, 0.051, 0.588, 0.236, 0.862, 0.022, 0.947, 0.145, 0.258, 0.595, 0.848, 0.633, 0.615, 0.014, 0.43, 0.529, 0.242, 0.011, 0.251, 0.524, 0.208, 0.011, 0.554, 0.271, 0.455, 0.021, 0.321]
global q = [0.059, 0.994, 0.873, 0.624, 0.882, 0.315, 0.983, 0.944, 0.258, 0.575, 0.954, 0.86, 0.852, 0.941, 0.917, 0.422, 0.987, 0.781, 0.773, 0.888, 0.367, 0.953, 0.617, 0.989, 0.363, 0.479, 0.927, 0.299, 0.887, 0.746, 0.922, 0.54, 0.718, 0.567, 0.488, 0.849, 0.654, 0.465, 0.888, 0.356, 0.992, 0.95, 0.703, 0.951, 0.479, 0.888, 0.86, 0.16, 0.653, 0.749, 0.962, 0.457, 0.998, 0.611, 0.782, 0.819, 0.889, 0.836, 0.771, 0.312, 0.963, 0.803, 0.73, 0.198, 0.929, 0.949, 0.834, 0.952, 0.981, 0.667, 0.91, 0.469, 0.714, 0.743, 0.927, 0.504, 0.467, 0.925, 0.345, 0.888, 0.619, 0.803, 0.692, 0.898, 0.271, 0.913, 0.94, 0.679, 0.972, 0.417, 0.936, 0.8, 0.356, 0.975, 0.971, 0.931, 0.314, 0.744, 0.88, 0.837, 0.494, 0.918, 0.829, 0.379, 0.921, 0.751, 0.941, 0.501, 0.589, 0.698, 0.947, 0.927, 0.915, 0.907, 0.893, 0.875, 0.985, 0.352, 0.933, 0.658, 0.779, 0.964, 0.756, 0.996, 0.666, 0.956, 0.228, 0.987, 0.33, 0.156, 0.518, 0.978, 0.321, 0.69, 0.163, 0.96, 0.866, 0.41, 0.735, 0.93, 0.56, 0.742, 0.861, 0.879, 0.812, 0.697, 0.502, 0.753, 0.802, 0.946, 0.724, 0.905, 0.438, 0.517, 0.946, 0.306, 0.987, 0.979, 0.82, 0.765, 0.999, 0.592, 0.965, 0.311, 0.822, 0.694, 0.83, 0.982, 0.742, 0.26, 0.666, 0.979, 0.866, 0.372, 0.468, 0.793, 0.613, 0.888, 0.823, 0.908, 0.986, 0.64, 0.781, 0.938, 0.394, 0.665, 0.115, 0.792, 0.611, 0.54, 0.913, 0.394, 0.554, 0.396, 0.437, 0.755, 0.902, 0.597, 0.915, 0.401, 0.881, 0.579, 0.988, 0.333, 0.949, 0.597, 0.558, 0.785, 0.97, 0.999, 0.855, 0.178, 0.697, 0.778, 0.533, 0.745, 0.471, 0.589, 0.251, 0.999, 0.805, 0.431, 0.806, 0.089, 0.831]
global origin = 1
global destination = 50