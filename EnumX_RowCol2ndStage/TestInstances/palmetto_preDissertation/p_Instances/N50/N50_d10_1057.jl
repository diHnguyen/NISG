global arcs = [1 14; 1 16; 1 32; 1 41; 1 46; 2 34; 2 37; 3 12; 3 23; 3 28; 3 31; 4 2; 4 5; 4 14; 4 20; 4 32; 4 38; 4 41; 4 43; 5 42; 6 8; 6 17; 6 25; 6 42; 7 9; 7 20; 7 29; 7 33; 7 44; 8 4; 8 11; 8 12; 9 11; 9 16; 9 18; 9 45; 9 49; 10 3; 10 4; 10 8; 10 21; 10 25; 10 50; 11 8; 11 9; 11 14; 11 17; 11 43; 11 47; 12 18; 12 25; 12 32; 12 44; 12 50; 13 15; 13 16; 13 18; 13 19; 13 28; 13 43; 13 46; 13 49; 14 12; 14 16; 14 17; 14 24; 14 28; 14 32; 14 37; 14 45; 14 47; 14 50; 15 2; 15 11; 15 26; 15 27; 15 33; 15 35; 16 2; 16 21; 16 40; 16 42; 16 50; 17 11; 17 28; 17 32; 17 33; 17 47; 18 2; 18 15; 18 22; 19 3; 19 11; 19 15; 19 17; 19 36; 19 37; 19 39; 20 10; 20 47; 21 5; 21 8; 21 19; 21 30; 21 36; 21 39; 21 40; 21 48; 22 7; 22 13; 22 16; 22 21; 22 45; 22 49; 22 50; 23 20; 23 34; 23 45; 23 48; 23 49; 24 16; 24 18; 24 20; 24 33; 25 4; 25 6; 25 12; 25 26; 25 33; 25 40; 25 41; 25 49; 26 16; 26 19; 26 23; 26 30; 26 36; 27 3; 27 13; 27 16; 27 22; 28 7; 28 9; 28 18; 28 33; 28 36; 28 45; 29 4; 29 43; 30 8; 30 17; 30 34; 30 35; 31 12; 31 40; 31 49; 31 50; 32 4; 32 5; 32 7; 32 11; 32 17; 32 18; 32 20; 32 28; 32 36; 33 10; 33 12; 33 23; 33 29; 33 32; 33 37; 33 43; 34 10; 34 16; 34 28; 34 33; 34 46; 34 47; 34 50; 35 7; 35 11; 35 12; 35 18; 35 27; 35 46; 35 48; 36 4; 36 14; 36 20; 36 29; 36 32; 36 35; 36 41; 36 42; 37 4; 37 19; 37 29; 37 41; 37 45; 37 49; 38 12; 38 33; 38 34; 39 14; 39 24; 39 25; 39 28; 39 32; 39 37; 39 50; 40 8; 40 22; 40 25; 40 32; 40 37; 40 44; 41 3; 41 36; 41 47; 41 48; 42 3; 42 6; 42 19; 42 20; 42 39; 43 50; 44 17; 44 19; 44 20; 44 33; 45 3; 45 6; 45 18; 45 33; 45 40; 45 46; 46 8; 46 26; 46 36; 46 39; 46 40; 47 38; 47 42; 48 14; 48 18; 48 24; 48 36; 48 46; 49 7; 49 8; 49 11; 49 27; 49 29; 49 32; 49 42]
global d_x = [7.0, 1.0, 4.0, 10.0, 4.0, 9.0, 1.0, 5.0, 8.0, 3.0, 2.0, 9.0, 3.0, 10.0, 2.0, 1.0, 6.0, 9.0, 10.0, 5.0, 4.0, 9.0, 7.0, 8.0, 3.0, 9.0, 1.0, 4.0, 10.0, 4.0, 9.0, 10.0, 2.0, 1.0, 10.0, 7.0, 1.0, 9.0, 7.0, 7.0, 3.0, 4.0, 8.0, 1.0, 4.0, 3.0, 8.0, 4.0, 1.0, 4.0, 10.0, 1.0, 10.0, 3.0, 1.0, 8.0, 2.0, 3.0, 4.0, 4.0, 9.0, 3.0, 4.0, 9.0, 4.0, 10.0, 6.0, 3.0, 5.0, 9.0, 8.0, 8.0, 4.0, 10.0, 8.0, 5.0, 6.0, 4.0, 5.0, 3.0, 9.0, 3.0, 8.0, 2.0, 3.0, 3.0, 7.0, 8.0, 4.0, 1.0, 5.0, 7.0, 2.0, 8.0, 6.0, 4.0, 8.0, 3.0, 3.0, 2.0, 4.0, 1.0, 3.0, 10.0, 5.0, 4.0, 3.0, 10.0, 10.0, 5.0, 3.0, 3.0, 8.0, 6.0, 4.0, 6.0, 1.0, 2.0, 1.0, 7.0, 4.0, 2.0, 5.0, 6.0, 9.0, 10.0, 4.0, 5.0, 2.0, 10.0, 5.0, 4.0, 9.0, 9.0, 1.0, 9.0, 7.0, 2.0, 8.0, 1.0, 5.0, 2.0, 2.0, 8.0, 1.0, 10.0, 3.0, 4.0, 10.0, 8.0, 3.0, 10.0, 9.0, 1.0, 7.0, 10.0, 10.0, 4.0, 9.0, 10.0, 8.0, 10.0, 6.0, 6.0, 5.0, 6.0, 9.0, 3.0, 9.0, 8.0, 10.0, 4.0, 1.0, 6.0, 5.0, 1.0, 8.0, 10.0, 10.0, 4.0, 9.0, 4.0, 10.0, 10.0, 1.0, 7.0, 3.0, 9.0, 1.0, 9.0, 8.0, 7.0, 8.0, 3.0, 9.0, 3.0, 1.0, 9.0, 1.0, 8.0, 5.0, 7.0, 3.0, 9.0, 2.0, 8.0, 5.0, 3.0, 5.0, 10.0, 1.0, 1.0, 9.0, 8.0, 2.0, 3.0, 8.0, 9.0, 7.0, 2.0, 9.0, 3.0, 8.0, 6.0, 3.0, 3.0, 1.0, 5.0, 4.0, 8.0, 7.0, 5.0, 6.0, 6.0, 2.0, 5.0, 8.0, 6.0, 6.0, 4.0, 2.0, 8.0, 2.0, 8.0, 6.0, 3.0, 2.0, 5.0, 8.0, 6.0, 8.0, 5.0, 3.0, 5.0, 1.0, 2.0]
global b_x = 5
global d_y = [8.0, 5.0, 8.0, 4.0, 7.0, 9.0, 7.0, 1.0, 2.0, 3.0, 4.0, 7.0, 9.0, 4.0, 8.0, 6.0, 10.0, 2.0, 2.0, 5.0, 9.0, 3.0, 9.0, 9.0, 1.0, 1.0, 5.0, 6.0, 10.0, 6.0, 8.0, 9.0, 7.0, 1.0, 6.0, 6.0, 3.0, 8.0, 6.0, 10.0, 9.0, 10.0, 2.0, 9.0, 2.0, 9.0, 5.0, 8.0, 4.0, 6.0, 8.0, 8.0, 8.0, 2.0, 5.0, 3.0, 1.0, 1.0, 10.0, 8.0, 7.0, 8.0, 10.0, 8.0, 6.0, 4.0, 8.0, 3.0, 1.0, 5.0, 2.0, 9.0, 3.0, 3.0, 8.0, 9.0, 7.0, 1.0, 3.0, 2.0, 4.0, 3.0, 4.0, 7.0, 6.0, 7.0, 4.0, 1.0, 5.0, 4.0, 2.0, 8.0, 5.0, 2.0, 3.0, 2.0, 10.0, 1.0, 9.0, 8.0, 2.0, 7.0, 10.0, 9.0, 7.0, 3.0, 7.0, 1.0, 10.0, 7.0, 4.0, 8.0, 5.0, 7.0, 10.0, 1.0, 1.0, 6.0, 5.0, 8.0, 9.0, 5.0, 7.0, 5.0, 5.0, 9.0, 10.0, 8.0, 3.0, 2.0, 6.0, 4.0, 9.0, 4.0, 4.0, 9.0, 6.0, 5.0, 2.0, 4.0, 2.0, 7.0, 4.0, 4.0, 6.0, 10.0, 4.0, 5.0, 6.0, 4.0, 6.0, 4.0, 2.0, 1.0, 5.0, 9.0, 2.0, 3.0, 10.0, 3.0, 3.0, 4.0, 9.0, 9.0, 2.0, 8.0, 1.0, 10.0, 6.0, 5.0, 4.0, 7.0, 2.0, 2.0, 1.0, 8.0, 10.0, 2.0, 10.0, 4.0, 3.0, 10.0, 7.0, 7.0, 8.0, 6.0, 10.0, 4.0, 1.0, 6.0, 3.0, 7.0, 6.0, 7.0, 3.0, 7.0, 4.0, 8.0, 4.0, 1.0, 4.0, 1.0, 10.0, 5.0, 4.0, 1.0, 7.0, 4.0, 2.0, 4.0, 5.0, 4.0, 5.0, 7.0, 1.0, 7.0, 6.0, 1.0, 10.0, 2.0, 5.0, 6.0, 4.0, 9.0, 3.0, 6.0, 8.0, 1.0, 4.0, 5.0, 8.0, 3.0, 6.0, 2.0, 7.0, 3.0, 8.0, 5.0, 9.0, 7.0, 9.0, 8.0, 5.0, 7.0, 1.0, 1.0, 4.0, 2.0, 8.0, 8.0, 10.0, 1.0, 6.0, 9.0, 6.0, 2.0]
global b_y = 10
global p = [0.744, 0.213, 0.272, 0.803, 0.204, 0.447, 0.162, 0.705, 0.76, 0.297, 0.355, 0.494, 0.634, 0.379, 0.46, 0.239, 0.614, 0.597, 0.315, 0.482, 0.81, 0.151, 0.821, 0.263, 0.928, 0.265, 0.045, 0.969, 0.37, 0.287, 0.23, 0.931, 0.919, 0.4, 0.632, 0.864, 0.863, 0.295, 0.744, 0.589, 0.428, 0.336, 0.135, 0.012, 0.346, 0.293, 0.304, 0.872, 0.177, 0.306, 0.181, 0.444, 0.688, 0.07, 0.679, 0.247, 0.623, 0.255, 0.197, 0.672, 0.657, 0.517, 0.146, 0.479, 0.485, 0.063, 0.31, 0.308, 0.185, 0.294, 0.571, 0.505, 0.088, 0.63, 0.871, 0.874, 0.808, 0.423, 0.095, 0.91, 0.972, 0.2, 0.586, 0.406, 0.206, 0.389, 0.71, 0.785, 0.554, 0.869, 0.966, 0.575, 0.344, 0.883, 0.907, 0.94, 0.741, 0.322, 0.153, 0.495, 0.336, 0.63, 0.071, 0.147, 0.695, 0.384, 0.89, 0.175, 0.812, 0.428, 0.963, 0.886, 0.738, 0.352, 0.449, 0.221, 0.812, 0.789, 0.304, 0.729, 0.791, 0.709, 0.908, 0.571, 0.309, 0.993, 0.059, 0.005, 0.054, 0.568, 0.449, 0.023, 0.214, 0.994, 0.012, 0.643, 0.231, 0.261, 0.262, 0.389, 0.059, 0.273, 0.254, 0.919, 0.919, 0.292, 0.775, 0.602, 0.55, 0.08, 0.562, 0.156, 0.962, 0.336, 0.52, 0.211, 0.898, 0.707, 0.619, 0.466, 0.515, 0.689, 0.275, 0.767, 0.366, 0.854, 0.655, 0.043, 0.144, 0.007, 0.612, 0.227, 0.82, 0.627, 0.32, 0.874, 0.023, 0.252, 0.878, 0.746, 0.861, 0.567, 0.229, 0.592, 0.539, 0.792, 0.86, 0.282, 0.761, 0.093, 0.609, 0.353, 0.489, 0.978, 0.897, 0.367, 0.802, 0.876, 0.454, 0.906, 0.923, 0.749, 0.161, 0.316, 0.32, 0.236, 0.703, 0.597, 0.341, 0.545, 0.035, 0.731, 0.935, 0.235, 0.318, 0.296, 0.71, 0.754, 0.948, 0.295, 0.059, 0.749, 0.528, 0.713, 0.413, 0.431, 0.547, 0.534, 0.514, 0.367, 0.602, 0.214, 0.556, 0.961, 0.736, 0.945, 0.87, 0.973, 0.968, 0.926, 0.913, 0.551, 0.301, 0.122, 0.796, 0.282, 0.286, 0.718, 0.264, 0.611, 0.649, 0.431, 0.402, 0.43, 0.035, 0.146]
global q = [0.982, 0.733, 0.557, 0.962, 0.651, 0.959, 0.314, 0.898, 0.881, 0.846, 0.704, 0.836, 0.812, 0.506, 0.704, 0.604, 0.894, 0.66, 0.454, 0.641, 0.819, 0.405, 0.925, 0.978, 0.941, 0.513, 0.39, 0.996, 0.852, 0.802, 0.92, 0.963, 0.94, 0.566, 0.808, 0.904, 0.91, 0.636, 0.961, 0.774, 0.544, 0.883, 0.587, 0.678, 0.944, 0.583, 0.343, 0.983, 0.462, 0.7, 0.184, 0.857, 0.758, 0.137, 0.982, 0.269, 0.92, 0.763, 0.533, 0.762, 0.969, 0.636, 0.175, 0.863, 0.846, 0.94, 0.654, 0.316, 0.897, 0.354, 0.774, 0.862, 0.395, 0.916, 0.945, 0.943, 0.982, 0.93, 0.804, 0.996, 0.995, 0.576, 0.834, 0.666, 0.683, 0.962, 0.91, 0.94, 0.789, 0.922, 0.991, 0.678, 0.583, 0.941, 0.959, 0.991, 0.903, 0.443, 0.522, 0.867, 0.502, 0.868, 0.153, 0.585, 0.964, 0.993, 0.951, 0.263, 0.862, 0.449, 0.998, 0.981, 0.997, 0.82, 0.513, 0.941, 0.89, 0.829, 0.369, 0.903, 0.932, 0.718, 0.988, 0.636, 0.67, 0.999, 0.185, 0.022, 0.545, 0.958, 0.882, 0.595, 0.424, 0.999, 0.098, 0.72, 0.555, 0.893, 0.512, 0.929, 0.573, 0.357, 0.317, 0.929, 0.952, 0.839, 0.971, 0.621, 0.74, 0.389, 0.762, 0.578, 0.999, 0.814, 0.797, 0.428, 0.981, 0.879, 0.841, 0.859, 0.716, 0.926, 0.697, 0.885, 0.753, 0.856, 0.9, 0.759, 0.824, 0.047, 0.828, 0.23, 0.869, 0.73, 0.892, 0.94, 0.625, 0.946, 0.924, 0.937, 0.93, 0.614, 0.763, 0.85, 0.835, 0.807, 0.916, 0.454, 0.856, 0.274, 0.624, 0.973, 0.942, 0.993, 0.941, 0.624, 0.954, 0.944, 0.693, 0.956, 0.938, 0.986, 0.517, 0.397, 0.941, 0.553, 0.855, 0.641, 0.404, 0.989, 0.423, 0.745, 0.94, 0.741, 0.451, 0.764, 0.882, 0.918, 0.992, 0.915, 0.926, 0.817, 0.797, 0.798, 0.732, 0.515, 0.836, 0.963, 0.69, 0.502, 0.885, 0.504, 0.752, 0.969, 0.749, 0.949, 0.986, 0.979, 0.98, 0.941, 0.932, 0.904, 0.947, 0.572, 0.903, 0.309, 0.508, 0.949, 0.39, 0.764, 0.672, 0.722, 0.966, 0.895, 0.634, 0.21]
global origin = 1
global destination = 50