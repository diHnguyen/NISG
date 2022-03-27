global arcs = [1 5; 1 11; 1 17; 1 18; 1 19; 1 23; 1 29; 1 45; 1 47; 2 5; 2 17; 2 21; 2 31; 2 42; 2 50; 3 15; 4 8; 4 13; 4 18; 4 19; 4 21; 4 43; 4 45; 5 2; 5 8; 5 35; 5 47; 6 23; 6 41; 7 6; 8 3; 8 4; 8 5; 8 7; 8 14; 8 17; 9 21; 9 27; 10 9; 10 13; 10 20; 11 10; 11 23; 11 32; 11 34; 12 21; 13 2; 13 5; 13 6; 13 8; 13 15; 13 24; 13 27; 13 30; 13 32; 13 33; 13 39; 13 42; 13 48; 14 4; 14 12; 14 21; 15 2; 15 40; 16 5; 16 24; 16 42; 16 46; 16 48; 17 21; 17 26; 18 11; 18 15; 18 21; 18 25; 18 28; 18 32; 18 41; 18 48; 19 9; 19 17; 19 30; 19 39; 20 3; 20 11; 20 19; 20 28; 20 33; 20 49; 21 2; 21 6; 21 16; 21 22; 21 25; 21 27; 21 28; 21 33; 21 34; 22 4; 22 7; 22 10; 22 18; 22 28; 22 33; 23 2; 23 5; 23 10; 23 14; 23 41; 23 48; 24 34; 25 18; 25 29; 25 32; 25 38; 25 39; 25 42; 25 43; 25 49; 25 50; 26 12; 26 13; 26 20; 27 12; 27 15; 27 21; 27 31; 27 39; 27 40; 27 46; 28 25; 28 34; 28 39; 28 45; 28 47; 29 25; 30 7; 30 25; 30 46; 31 21; 31 27; 31 33; 31 37; 32 6; 32 19; 32 47; 33 14; 33 19; 33 26; 33 46; 33 50; 34 18; 34 30; 34 43; 34 47; 34 48; 35 9; 35 13; 35 19; 36 6; 36 29; 36 41; 36 43; 36 47; 37 4; 37 9; 37 12; 37 48; 38 8; 38 17; 38 23; 38 46; 38 49; 39 10; 39 11; 39 14; 39 38; 39 50; 40 6; 40 8; 40 21; 40 39; 40 46; 40 48; 40 49; 41 5; 41 36; 41 42; 42 14; 42 24; 42 32; 42 43; 42 44; 42 46; 43 9; 43 19; 44 3; 44 24; 44 34; 44 35; 44 41; 44 48; 45 14; 45 26; 45 29; 45 34; 45 41; 46 25; 46 36; 46 41; 46 48; 46 50; 47 3; 47 13; 47 16; 47 17; 48 9; 48 10; 48 11; 48 17; 48 27; 49 2; 49 13; 49 17; 49 22; 49 23; 49 24; 49 29]
global d_x = [6.0, 7.0, 2.0, 1.0, 1.0, 9.0, 2.0, 7.0, 5.0, 4.0, 7.0, 10.0, 2.0, 3.0, 1.0, 6.0, 4.0, 3.0, 9.0, 4.0, 3.0, 9.0, 5.0, 1.0, 4.0, 8.0, 1.0, 2.0, 4.0, 10.0, 9.0, 6.0, 6.0, 8.0, 10.0, 5.0, 3.0, 8.0, 2.0, 7.0, 4.0, 2.0, 2.0, 4.0, 8.0, 4.0, 7.0, 6.0, 2.0, 4.0, 6.0, 8.0, 6.0, 1.0, 6.0, 1.0, 2.0, 9.0, 1.0, 7.0, 4.0, 2.0, 1.0, 3.0, 6.0, 1.0, 10.0, 9.0, 6.0, 6.0, 4.0, 6.0, 4.0, 3.0, 6.0, 6.0, 5.0, 2.0, 2.0, 3.0, 2.0, 1.0, 8.0, 4.0, 1.0, 1.0, 9.0, 1.0, 10.0, 1.0, 4.0, 2.0, 2.0, 3.0, 4.0, 3.0, 4.0, 5.0, 7.0, 2.0, 4.0, 9.0, 2.0, 7.0, 2.0, 9.0, 10.0, 10.0, 1.0, 10.0, 10.0, 7.0, 8.0, 8.0, 8.0, 1.0, 9.0, 8.0, 4.0, 3.0, 3.0, 1.0, 9.0, 6.0, 2.0, 3.0, 4.0, 10.0, 4.0, 2.0, 7.0, 8.0, 4.0, 1.0, 9.0, 8.0, 9.0, 9.0, 9.0, 1.0, 5.0, 3.0, 3.0, 1.0, 10.0, 4.0, 5.0, 2.0, 7.0, 8.0, 10.0, 2.0, 3.0, 3.0, 10.0, 4.0, 5.0, 5.0, 4.0, 9.0, 6.0, 3.0, 10.0, 4.0, 6.0, 9.0, 10.0, 5.0, 9.0, 4.0, 3.0, 7.0, 7.0, 1.0, 9.0, 10.0, 4.0, 4.0, 5.0, 4.0, 3.0, 5.0, 2.0, 9.0, 10.0, 9.0, 4.0, 2.0, 7.0, 10.0, 5.0, 6.0, 4.0, 9.0, 7.0, 6.0, 1.0, 9.0, 7.0, 3.0, 6.0, 7.0, 3.0, 8.0, 2.0, 3.0, 2.0, 3.0, 6.0, 3.0, 6.0, 7.0, 6.0, 6.0, 7.0, 2.0, 2.0, 2.0, 10.0, 8.0, 2.0, 9.0, 3.0, 1.0, 3.0, 2.0, 3.0, 5.0]
global b_x = 5
global d_y = [6.0, 5.0, 9.0, 1.0, 4.0, 4.0, 10.0, 1.0, 7.0, 8.0, 7.0, 3.0, 1.0, 5.0, 5.0, 5.0, 9.0, 5.0, 9.0, 3.0, 2.0, 4.0, 2.0, 10.0, 7.0, 8.0, 4.0, 10.0, 5.0, 7.0, 6.0, 9.0, 3.0, 2.0, 2.0, 4.0, 9.0, 1.0, 6.0, 3.0, 6.0, 9.0, 8.0, 4.0, 7.0, 8.0, 7.0, 3.0, 4.0, 6.0, 7.0, 3.0, 4.0, 5.0, 6.0, 6.0, 3.0, 3.0, 10.0, 8.0, 5.0, 2.0, 7.0, 3.0, 4.0, 9.0, 9.0, 2.0, 2.0, 7.0, 7.0, 2.0, 9.0, 7.0, 6.0, 2.0, 6.0, 4.0, 1.0, 2.0, 1.0, 6.0, 6.0, 3.0, 7.0, 9.0, 10.0, 6.0, 6.0, 2.0, 7.0, 7.0, 3.0, 3.0, 3.0, 8.0, 8.0, 5.0, 9.0, 2.0, 8.0, 7.0, 2.0, 1.0, 5.0, 3.0, 2.0, 1.0, 5.0, 9.0, 10.0, 4.0, 1.0, 9.0, 1.0, 2.0, 9.0, 10.0, 4.0, 3.0, 5.0, 2.0, 10.0, 2.0, 9.0, 6.0, 7.0, 1.0, 5.0, 5.0, 5.0, 3.0, 4.0, 8.0, 2.0, 10.0, 1.0, 2.0, 2.0, 9.0, 8.0, 5.0, 5.0, 7.0, 7.0, 8.0, 8.0, 8.0, 8.0, 10.0, 7.0, 10.0, 4.0, 5.0, 4.0, 7.0, 9.0, 7.0, 9.0, 8.0, 6.0, 7.0, 10.0, 10.0, 9.0, 2.0, 7.0, 7.0, 9.0, 1.0, 4.0, 2.0, 4.0, 4.0, 2.0, 7.0, 5.0, 1.0, 2.0, 5.0, 10.0, 3.0, 3.0, 8.0, 5.0, 10.0, 8.0, 6.0, 8.0, 8.0, 3.0, 2.0, 2.0, 10.0, 1.0, 9.0, 5.0, 4.0, 9.0, 5.0, 6.0, 2.0, 4.0, 5.0, 7.0, 1.0, 5.0, 5.0, 5.0, 10.0, 7.0, 4.0, 8.0, 8.0, 7.0, 8.0, 1.0, 2.0, 9.0, 6.0, 8.0, 5.0, 7.0, 3.0, 5.0, 10.0, 7.0, 7.0]
global b_y = 10
global p = [0.59, 0.389, 0.394, 0.265, 0.297, 0.08, 0.534, 0.771, 0.599, 0.248, 0.515, 0.89, 0.614, 0.744, 0.693, 0.905, 0.414, 0.56, 0.742, 0.605, 0.676, 0.255, 0.077, 0.546, 0.359, 0.349, 0.225, 0.961, 0.045, 0.466, 0.621, 0.321, 0.273, 0.976, 0.67, 0.609, 0.836, 0.176, 0.536, 0.526, 0.192, 0.312, 0.855, 0.901, 0.339, 0.991, 0.882, 0.779, 0.187, 0.771, 0.68, 0.65, 0.561, 0.902, 0.013, 0.629, 0.291, 0.747, 0.288, 0.647, 0.184, 0.976, 0.689, 0.652, 0.619, 0.517, 0.05, 0.171, 0.847, 0.223, 0.284, 0.431, 0.242, 0.435, 0.671, 0.13, 0.941, 0.925, 0.33, 0.247, 0.178, 0.239, 0.811, 0.166, 0.498, 0.547, 0.292, 0.707, 0.484, 0.381, 0.553, 0.091, 0.15, 0.052, 0.146, 0.431, 0.462, 0.297, 0.279, 0.476, 0.27, 0.038, 0.039, 0.838, 0.851, 0.931, 0.078, 0.922, 0.721, 0.269, 0.144, 0.317, 0.847, 0.65, 0.634, 0.621, 0.378, 0.969, 0.437, 0.801, 0.472, 0.774, 0.885, 0.851, 0.077, 0.037, 0.933, 0.763, 0.329, 0.335, 0.481, 0.552, 0.743, 0.663, 0.335, 0.159, 0.059, 0.86, 0.972, 0.074, 0.118, 0.637, 0.572, 0.196, 0.116, 0.746, 0.821, 0.176, 0.235, 0.363, 0.379, 0.036, 0.319, 0.911, 0.761, 0.806, 0.222, 0.797, 0.459, 0.841, 0.063, 0.896, 0.048, 0.828, 0.671, 0.357, 0.913, 0.229, 0.706, 0.72, 0.202, 0.045, 0.212, 0.218, 0.028, 0.422, 0.894, 0.107, 0.738, 0.654, 0.402, 0.75, 0.112, 0.974, 0.483, 0.778, 0.627, 0.338, 0.388, 0.347, 0.181, 0.344, 0.363, 0.51, 0.461, 0.27, 0.681, 0.705, 0.668, 0.801, 0.23, 0.32, 0.862, 0.037, 0.431, 0.702, 0.809, 0.761, 0.713, 0.14, 0.609, 0.407, 0.566, 0.933, 0.336, 0.135, 0.72, 0.121, 0.515, 0.709, 0.599, 0.157, 0.592, 0.85, 0.663, 0.677, 0.338, 0.903]
global q = [0.654, 0.944, 0.758, 0.288, 0.763, 0.222, 0.852, 0.962, 0.918, 0.997, 0.671, 0.908, 0.988, 0.795, 0.814, 0.981, 0.693, 0.867, 0.746, 0.886, 0.866, 0.749, 0.83, 0.947, 0.509, 0.615, 0.419, 0.967, 0.535, 0.523, 0.964, 0.852, 0.907, 0.999, 0.952, 0.712, 0.983, 0.255, 0.778, 0.86, 0.482, 0.537, 0.905, 0.915, 0.785, 0.992, 0.944, 0.912, 0.188, 0.941, 0.953, 0.816, 0.623, 0.959, 0.57, 0.801, 0.984, 0.78, 0.934, 0.751, 0.288, 0.981, 0.969, 0.765, 0.711, 0.755, 0.862, 0.958, 0.981, 0.586, 0.723, 0.556, 0.966, 0.6, 0.733, 0.851, 0.969, 0.976, 0.807, 0.526, 0.79, 0.464, 0.837, 0.193, 0.984, 0.556, 0.744, 0.871, 0.679, 0.583, 0.703, 0.469, 0.988, 0.234, 0.757, 0.998, 0.99, 0.677, 0.639, 0.689, 0.3, 0.81, 0.726, 0.946, 0.904, 0.968, 0.238, 0.998, 0.921, 0.577, 0.437, 0.457, 0.851, 0.653, 0.828, 0.995, 0.553, 0.992, 0.483, 0.867, 0.8, 0.863, 0.92, 0.948, 0.308, 0.819, 0.954, 0.902, 0.737, 0.86, 0.683, 0.736, 0.92, 0.842, 0.584, 0.201, 0.188, 0.931, 0.985, 0.519, 0.434, 0.652, 0.697, 0.956, 0.266, 0.809, 0.883, 0.547, 0.539, 0.613, 0.627, 0.418, 0.734, 0.929, 0.89, 0.858, 0.517, 0.865, 0.678, 0.88, 0.981, 0.937, 0.17, 0.945, 0.725, 0.826, 0.925, 0.546, 0.857, 0.851, 0.583, 0.082, 0.325, 0.588, 0.043, 0.856, 0.943, 0.109, 0.932, 0.827, 0.954, 0.805, 0.846, 0.981, 0.905, 0.843, 0.742, 0.772, 0.524, 0.95, 0.553, 0.912, 0.502, 0.977, 0.516, 0.454, 0.922, 0.706, 0.744, 0.97, 0.388, 0.543, 0.883, 0.25, 0.438, 0.747, 0.994, 0.804, 0.725, 0.304, 0.729, 0.693, 0.71, 0.993, 0.639, 0.483, 0.926, 0.641, 0.758, 0.861, 0.697, 0.733, 0.718, 0.966, 0.739, 0.84, 0.52, 0.918]
global origin = 1
global destination = 50