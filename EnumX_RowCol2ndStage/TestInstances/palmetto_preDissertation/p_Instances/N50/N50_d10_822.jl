global arcs = [1 7; 1 12; 1 28; 1 44; 1 50; 2 9; 2 27; 2 31; 2 40; 2 42; 2 45; 2 48; 3 14; 3 18; 3 19; 3 27; 4 8; 4 19; 4 25; 4 30; 4 44; 5 4; 5 17; 5 18; 5 30; 5 33; 5 42; 5 45; 5 49; 6 13; 6 21; 6 25; 6 40; 6 46; 7 10; 7 24; 7 27; 7 29; 7 38; 7 42; 7 44; 7 47; 8 29; 8 30; 8 34; 8 44; 8 45; 9 6; 9 8; 9 24; 9 33; 9 47; 10 6; 10 8; 10 12; 10 13; 10 43; 11 2; 11 18; 11 28; 11 50; 12 9; 12 14; 12 18; 12 39; 12 40; 12 45; 13 15; 13 31; 13 35; 14 2; 14 28; 14 34; 14 41; 15 3; 15 38; 15 46; 16 23; 16 33; 17 2; 17 7; 17 15; 17 33; 17 41; 17 46; 17 50; 18 2; 18 16; 18 37; 18 39; 18 49; 18 50; 19 5; 19 24; 19 29; 19 36; 19 46; 19 50; 20 5; 20 14; 20 44; 21 29; 21 39; 21 40; 21 48; 22 3; 22 7; 22 8; 22 46; 23 2; 23 9; 23 35; 24 21; 24 23; 24 26; 24 27; 24 31; 24 33; 24 43; 25 20; 25 21; 25 26; 25 31; 26 8; 26 9; 26 17; 26 20; 26 22; 26 24; 26 28; 26 43; 26 47; 27 17; 27 18; 27 24; 27 36; 27 43; 28 7; 28 9; 28 13; 28 27; 28 31; 28 33; 29 3; 29 10; 29 21; 29 25; 29 30; 29 44; 29 45; 30 13; 30 44; 31 2; 31 11; 31 25; 31 37; 32 9; 32 14; 32 29; 32 39; 33 14; 33 17; 33 29; 34 4; 34 22; 34 44; 34 49; 35 4; 35 21; 35 24; 35 25; 35 27; 35 31; 35 38; 35 41; 36 20; 36 21; 37 5; 37 25; 37 28; 38 7; 38 13; 38 17; 38 37; 39 9; 39 28; 39 38; 39 43; 39 46; 40 12; 40 30; 40 50; 41 36; 42 14; 42 15; 42 23; 42 29; 43 6; 43 14; 43 20; 43 29; 43 32; 44 16; 44 29; 44 42; 45 5; 46 12; 46 22; 46 42; 47 4; 47 21; 47 31; 47 35; 47 39; 47 44; 47 45; 48 9; 48 26; 48 33; 48 38; 49 13; 49 23; 49 27; 49 32; 49 34; 49 39; 49 46]
global d_x = [6.0, 9.0, 5.0, 4.0, 1.0, 5.0, 9.0, 2.0, 3.0, 5.0, 5.0, 2.0, 3.0, 9.0, 9.0, 8.0, 6.0, 8.0, 1.0, 2.0, 4.0, 7.0, 1.0, 9.0, 3.0, 8.0, 2.0, 8.0, 4.0, 5.0, 2.0, 9.0, 8.0, 1.0, 9.0, 2.0, 1.0, 9.0, 9.0, 10.0, 5.0, 2.0, 4.0, 9.0, 1.0, 2.0, 3.0, 7.0, 2.0, 8.0, 2.0, 4.0, 9.0, 3.0, 1.0, 1.0, 1.0, 1.0, 7.0, 9.0, 5.0, 4.0, 7.0, 7.0, 1.0, 1.0, 1.0, 2.0, 3.0, 9.0, 9.0, 10.0, 5.0, 8.0, 1.0, 7.0, 8.0, 10.0, 3.0, 8.0, 6.0, 8.0, 2.0, 8.0, 8.0, 3.0, 8.0, 5.0, 2.0, 1.0, 8.0, 2.0, 10.0, 4.0, 4.0, 5.0, 8.0, 7.0, 6.0, 8.0, 7.0, 5.0, 4.0, 1.0, 10.0, 7.0, 2.0, 6.0, 1.0, 1.0, 3.0, 6.0, 7.0, 2.0, 1.0, 3.0, 1.0, 9.0, 5.0, 9.0, 6.0, 10.0, 8.0, 2.0, 8.0, 1.0, 2.0, 1.0, 7.0, 4.0, 9.0, 6.0, 3.0, 7.0, 8.0, 2.0, 10.0, 2.0, 2.0, 3.0, 8.0, 10.0, 2.0, 7.0, 7.0, 8.0, 3.0, 5.0, 1.0, 8.0, 4.0, 5.0, 1.0, 3.0, 6.0, 1.0, 1.0, 9.0, 6.0, 5.0, 7.0, 1.0, 4.0, 1.0, 4.0, 2.0, 6.0, 4.0, 3.0, 10.0, 1.0, 5.0, 2.0, 9.0, 5.0, 10.0, 1.0, 7.0, 5.0, 10.0, 10.0, 5.0, 1.0, 9.0, 3.0, 4.0, 2.0, 6.0, 9.0, 7.0, 2.0, 8.0, 3.0, 7.0, 4.0, 2.0, 3.0, 6.0, 5.0, 8.0, 3.0, 3.0, 9.0, 2.0, 4.0, 2.0, 5.0, 1.0, 5.0, 6.0, 2.0, 2.0, 8.0, 10.0, 9.0, 6.0, 2.0, 10.0, 3.0, 8.0, 9.0, 8.0, 8.0, 4.0, 2.0, 5.0, 6.0]
global b_x = 5
global d_y = [3.0, 9.0, 10.0, 6.0, 5.0, 2.0, 1.0, 5.0, 2.0, 2.0, 6.0, 3.0, 5.0, 7.0, 9.0, 2.0, 2.0, 3.0, 5.0, 1.0, 4.0, 2.0, 1.0, 6.0, 5.0, 1.0, 1.0, 1.0, 9.0, 2.0, 3.0, 5.0, 1.0, 8.0, 4.0, 5.0, 7.0, 1.0, 8.0, 7.0, 4.0, 4.0, 8.0, 10.0, 5.0, 5.0, 5.0, 8.0, 10.0, 1.0, 1.0, 7.0, 7.0, 7.0, 8.0, 5.0, 5.0, 4.0, 8.0, 6.0, 8.0, 4.0, 4.0, 4.0, 4.0, 6.0, 8.0, 4.0, 1.0, 1.0, 4.0, 4.0, 7.0, 5.0, 9.0, 3.0, 7.0, 1.0, 4.0, 3.0, 10.0, 2.0, 4.0, 10.0, 8.0, 5.0, 5.0, 2.0, 5.0, 5.0, 6.0, 6.0, 5.0, 9.0, 8.0, 2.0, 8.0, 1.0, 7.0, 5.0, 9.0, 3.0, 8.0, 5.0, 9.0, 10.0, 7.0, 1.0, 7.0, 5.0, 9.0, 7.0, 8.0, 10.0, 1.0, 1.0, 6.0, 4.0, 8.0, 7.0, 7.0, 4.0, 9.0, 3.0, 8.0, 8.0, 2.0, 1.0, 9.0, 9.0, 7.0, 3.0, 9.0, 9.0, 2.0, 8.0, 10.0, 6.0, 1.0, 8.0, 5.0, 2.0, 3.0, 8.0, 1.0, 8.0, 8.0, 8.0, 2.0, 7.0, 5.0, 4.0, 7.0, 7.0, 8.0, 1.0, 6.0, 2.0, 2.0, 7.0, 4.0, 6.0, 7.0, 7.0, 9.0, 10.0, 7.0, 5.0, 10.0, 9.0, 3.0, 6.0, 3.0, 7.0, 10.0, 10.0, 4.0, 6.0, 4.0, 5.0, 10.0, 8.0, 9.0, 4.0, 10.0, 8.0, 5.0, 3.0, 10.0, 9.0, 5.0, 10.0, 1.0, 3.0, 5.0, 1.0, 10.0, 2.0, 3.0, 3.0, 10.0, 4.0, 8.0, 7.0, 6.0, 2.0, 5.0, 3.0, 8.0, 3.0, 7.0, 10.0, 6.0, 8.0, 9.0, 3.0, 7.0, 4.0, 2.0, 8.0, 8.0, 3.0, 9.0, 8.0, 2.0, 4.0, 4.0]
global b_y = 10
global p = [0.07, 0.743, 0.648, 0.298, 0.375, 0.818, 0.511, 0.77, 0.354, 0.762, 0.185, 0.686, 0.349, 0.631, 0.943, 0.452, 0.785, 0.866, 0.23, 0.813, 0.087, 0.707, 0.8, 0.028, 0.176, 0.019, 0.62, 0.312, 0.949, 0.106, 0.723, 0.437, 0.47, 0.276, 0.865, 0.408, 0.97, 0.604, 0.381, 0.799, 0.179, 0.707, 0.494, 0.321, 0.202, 0.224, 0.643, 0.463, 0.167, 0.684, 0.06, 0.293, 0.876, 0.092, 0.106, 0.814, 0.794, 0.811, 0.086, 0.955, 0.25, 0.462, 0.161, 0.163, 0.61, 0.237, 0.504, 0.52, 0.945, 0.568, 0.693, 0.938, 0.688, 0.657, 0.22, 0.281, 0.004, 0.948, 0.323, 0.2, 0.705, 0.106, 0.15, 0.049, 0.402, 0.505, 0.645, 0.289, 0.99, 0.692, 0.317, 0.078, 0.706, 0.341, 0.013, 0.896, 0.818, 0.45, 0.583, 0.293, 0.317, 0.648, 0.182, 0.619, 0.475, 0.323, 0.053, 0.984, 0.428, 0.533, 0.894, 0.601, 0.718, 0.44, 0.713, 0.48, 0.4, 0.265, 0.231, 0.911, 0.36, 0.602, 0.784, 0.838, 0.776, 0.757, 0.941, 0.478, 0.875, 0.553, 0.153, 0.955, 0.358, 0.573, 0.136, 0.363, 0.429, 0.58, 0.988, 0.822, 0.639, 0.914, 0.723, 0.707, 0.289, 0.99, 0.104, 0.332, 0.511, 0.029, 0.681, 0.048, 0.117, 0.35, 0.713, 0.231, 0.476, 0.659, 0.017, 0.39, 0.472, 0.677, 0.393, 0.494, 0.712, 0.359, 0.966, 0.013, 0.981, 0.158, 0.164, 0.728, 0.75, 0.247, 0.291, 0.179, 0.248, 0.544, 0.301, 0.323, 0.474, 0.814, 0.72, 0.405, 0.878, 0.176, 0.832, 0.511, 0.756, 0.879, 0.892, 0.038, 0.056, 0.462, 0.515, 0.688, 0.363, 0.04, 0.261, 0.788, 0.94, 0.99, 0.437, 0.188, 0.031, 0.631, 0.767, 0.067, 0.834, 0.409, 0.493, 0.317, 0.515, 0.906, 0.741, 0.507, 0.483, 0.247, 0.761, 0.728, 0.384, 0.553, 0.07, 0.109, 0.009, 0.413, 0.255]
global q = [0.816, 0.81, 0.938, 0.841, 0.821, 0.83, 0.579, 0.885, 0.871, 0.988, 0.363, 0.828, 0.63, 0.846, 0.988, 0.821, 0.889, 0.884, 0.641, 0.966, 0.213, 0.808, 0.807, 0.6, 0.967, 0.758, 0.64, 0.414, 0.988, 0.959, 0.931, 0.816, 0.897, 0.377, 0.94, 0.444, 0.977, 0.839, 0.713, 0.816, 0.69, 0.749, 0.765, 0.83, 0.421, 0.398, 0.794, 0.551, 0.283, 0.688, 0.791, 0.631, 0.904, 0.598, 0.77, 0.969, 0.948, 0.947, 0.863, 0.955, 0.555, 0.748, 0.546, 0.561, 0.805, 0.716, 0.853, 0.621, 0.957, 0.882, 0.831, 0.959, 0.731, 0.953, 0.416, 0.805, 0.299, 0.961, 0.665, 0.845, 0.841, 0.961, 0.427, 0.734, 0.928, 0.726, 0.989, 0.476, 0.991, 0.953, 0.455, 0.444, 0.81, 0.846, 0.363, 0.968, 0.986, 0.547, 0.933, 0.551, 0.65, 0.725, 0.972, 0.868, 0.748, 0.65, 0.252, 0.985, 0.764, 0.556, 0.894, 0.8, 0.943, 0.991, 0.894, 0.902, 0.553, 0.377, 0.821, 0.99, 0.518, 0.725, 0.798, 0.951, 0.988, 0.996, 0.996, 0.954, 0.937, 0.576, 0.454, 0.958, 0.971, 0.875, 0.538, 0.786, 0.738, 0.695, 0.992, 0.823, 0.996, 0.954, 0.837, 0.875, 0.459, 0.996, 0.995, 0.405, 0.682, 0.134, 0.715, 0.406, 0.927, 0.516, 0.784, 0.476, 0.753, 0.739, 0.892, 0.768, 0.583, 0.696, 0.409, 0.95, 0.856, 0.567, 0.974, 0.117, 0.987, 0.399, 0.665, 0.818, 0.99, 0.971, 0.679, 0.683, 0.423, 0.786, 0.917, 0.669, 0.816, 0.86, 0.811, 0.805, 0.908, 0.747, 0.935, 0.71, 0.903, 0.966, 0.934, 0.553, 0.906, 0.687, 0.906, 0.884, 0.991, 0.642, 0.673, 0.99, 0.969, 0.996, 0.998, 0.712, 0.195, 0.949, 0.788, 0.133, 0.943, 0.731, 0.973, 0.681, 0.588, 0.99, 0.959, 0.602, 0.608, 0.332, 0.904, 0.942, 0.499, 0.83, 0.44, 0.664, 0.248, 0.733, 0.307]
global origin = 1
global destination = 50