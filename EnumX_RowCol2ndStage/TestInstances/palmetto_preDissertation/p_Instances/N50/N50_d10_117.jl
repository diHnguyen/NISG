global arcs = [1 9; 1 28; 1 31; 1 43; 2 24; 2 30; 2 33; 2 35; 2 46; 2 49; 3 11; 3 29; 4 12; 4 19; 4 20; 4 34; 4 40; 4 47; 5 15; 6 23; 6 28; 6 43; 6 48; 7 3; 7 8; 7 22; 7 34; 7 35; 7 39; 8 26; 8 38; 8 40; 8 42; 8 48; 9 14; 9 16; 9 39; 9 48; 10 33; 10 43; 10 45; 11 6; 11 33; 11 35; 11 38; 11 49; 12 7; 12 11; 12 16; 12 28; 12 30; 12 34; 12 48; 13 17; 13 26; 13 29; 13 30; 13 45; 13 48; 14 24; 14 37; 14 43; 15 2; 15 7; 15 12; 15 20; 15 22; 15 30; 15 38; 15 43; 15 45; 16 27; 16 38; 16 46; 17 9; 17 18; 17 27; 17 30; 17 32; 17 38; 17 43; 17 46; 17 48; 17 50; 18 16; 18 21; 18 41; 19 7; 19 13; 19 21; 19 26; 19 27; 20 2; 20 7; 20 8; 20 30; 20 38; 20 43; 20 49; 21 24; 21 37; 21 42; 21 44; 21 48; 22 2; 22 13; 22 20; 22 21; 22 23; 22 38; 23 4; 23 6; 23 14; 23 16; 23 22; 23 24; 23 35; 23 49; 24 10; 24 12; 24 25; 25 8; 25 13; 25 20; 25 48; 26 13; 26 20; 26 36; 27 9; 28 3; 28 4; 28 7; 28 15; 28 32; 28 42; 29 9; 29 20; 29 25; 29 40; 29 41; 29 43; 29 44; 29 46; 30 31; 30 37; 30 39; 30 45; 30 48; 31 9; 31 33; 31 38; 31 50; 32 10; 32 17; 32 23; 32 39; 32 48; 32 49; 33 4; 33 22; 33 26; 33 35; 33 37; 33 44; 33 47; 34 24; 34 28; 34 43; 35 16; 35 28; 35 40; 35 48; 36 4; 36 20; 36 23; 36 29; 36 38; 36 45; 37 18; 37 19; 37 41; 38 17; 38 25; 38 26; 38 35; 38 41; 38 45; 39 4; 39 6; 39 10; 39 11; 39 21; 40 3; 40 5; 40 16; 40 17; 40 24; 41 9; 41 18; 41 21; 41 31; 42 12; 42 35; 42 37; 42 38; 43 2; 43 3; 43 5; 43 7; 43 10; 43 17; 43 37; 43 44; 44 18; 44 24; 44 29; 44 40; 44 50; 45 9; 45 26; 45 44; 45 50; 46 32; 46 50; 47 27; 47 32; 47 34; 47 35; 48 6; 48 10; 48 16; 48 18; 48 19; 48 26; 48 30; 49 6; 49 12; 49 35]
global d_x = [1.0, 3.0, 7.0, 10.0, 5.0, 8.0, 7.0, 1.0, 5.0, 5.0, 1.0, 4.0, 10.0, 4.0, 4.0, 9.0, 9.0, 3.0, 7.0, 6.0, 5.0, 3.0, 5.0, 2.0, 3.0, 8.0, 9.0, 4.0, 10.0, 3.0, 8.0, 4.0, 6.0, 4.0, 3.0, 2.0, 4.0, 5.0, 1.0, 6.0, 3.0, 3.0, 4.0, 4.0, 5.0, 9.0, 10.0, 4.0, 2.0, 7.0, 9.0, 8.0, 9.0, 5.0, 8.0, 2.0, 6.0, 4.0, 8.0, 2.0, 7.0, 7.0, 2.0, 9.0, 8.0, 5.0, 10.0, 2.0, 2.0, 10.0, 1.0, 4.0, 2.0, 5.0, 7.0, 5.0, 2.0, 10.0, 6.0, 4.0, 5.0, 3.0, 1.0, 4.0, 9.0, 6.0, 10.0, 7.0, 4.0, 5.0, 6.0, 6.0, 8.0, 7.0, 4.0, 8.0, 5.0, 4.0, 4.0, 7.0, 3.0, 4.0, 4.0, 1.0, 2.0, 2.0, 6.0, 5.0, 1.0, 8.0, 2.0, 10.0, 3.0, 2.0, 8.0, 5.0, 9.0, 5.0, 4.0, 4.0, 6.0, 6.0, 3.0, 9.0, 5.0, 5.0, 2.0, 7.0, 2.0, 9.0, 5.0, 5.0, 10.0, 4.0, 1.0, 7.0, 3.0, 8.0, 4.0, 3.0, 10.0, 2.0, 2.0, 6.0, 5.0, 10.0, 6.0, 9.0, 1.0, 10.0, 7.0, 2.0, 2.0, 3.0, 3.0, 8.0, 5.0, 8.0, 3.0, 10.0, 1.0, 8.0, 3.0, 8.0, 10.0, 5.0, 9.0, 1.0, 2.0, 2.0, 1.0, 10.0, 4.0, 6.0, 6.0, 2.0, 6.0, 9.0, 2.0, 6.0, 3.0, 8.0, 5.0, 1.0, 7.0, 7.0, 7.0, 6.0, 2.0, 10.0, 4.0, 6.0, 8.0, 10.0, 5.0, 7.0, 10.0, 1.0, 5.0, 4.0, 4.0, 8.0, 2.0, 7.0, 8.0, 3.0, 7.0, 4.0, 3.0, 4.0, 1.0, 3.0, 1.0, 3.0, 3.0, 6.0, 9.0, 7.0, 2.0, 9.0, 8.0, 3.0, 10.0, 1.0, 7.0, 5.0, 1.0, 3.0, 2.0, 9.0, 10.0, 7.0, 8.0, 1.0, 2.0, 8.0, 10.0, 4.0]
global b_x = 5
global d_y = [4.0, 10.0, 6.0, 6.0, 6.0, 5.0, 5.0, 3.0, 8.0, 6.0, 10.0, 3.0, 10.0, 10.0, 1.0, 8.0, 2.0, 7.0, 6.0, 10.0, 8.0, 9.0, 10.0, 3.0, 6.0, 4.0, 3.0, 7.0, 4.0, 8.0, 6.0, 1.0, 7.0, 2.0, 2.0, 5.0, 10.0, 3.0, 5.0, 8.0, 7.0, 9.0, 7.0, 7.0, 3.0, 2.0, 1.0, 10.0, 7.0, 5.0, 9.0, 9.0, 6.0, 4.0, 5.0, 9.0, 10.0, 2.0, 10.0, 6.0, 6.0, 2.0, 9.0, 8.0, 3.0, 7.0, 1.0, 1.0, 5.0, 1.0, 2.0, 10.0, 9.0, 8.0, 6.0, 6.0, 4.0, 10.0, 1.0, 6.0, 1.0, 3.0, 7.0, 2.0, 6.0, 6.0, 10.0, 5.0, 10.0, 2.0, 6.0, 6.0, 4.0, 10.0, 5.0, 6.0, 10.0, 9.0, 3.0, 5.0, 10.0, 7.0, 9.0, 9.0, 2.0, 2.0, 6.0, 4.0, 5.0, 10.0, 6.0, 4.0, 2.0, 6.0, 6.0, 1.0, 2.0, 9.0, 2.0, 4.0, 9.0, 2.0, 2.0, 8.0, 3.0, 1.0, 1.0, 3.0, 2.0, 4.0, 2.0, 8.0, 4.0, 4.0, 1.0, 3.0, 4.0, 6.0, 1.0, 4.0, 10.0, 10.0, 2.0, 8.0, 1.0, 10.0, 1.0, 8.0, 8.0, 10.0, 7.0, 4.0, 8.0, 4.0, 6.0, 6.0, 8.0, 6.0, 10.0, 8.0, 7.0, 8.0, 9.0, 1.0, 3.0, 9.0, 3.0, 2.0, 4.0, 4.0, 8.0, 6.0, 4.0, 1.0, 7.0, 1.0, 10.0, 4.0, 3.0, 9.0, 8.0, 8.0, 1.0, 6.0, 3.0, 4.0, 10.0, 3.0, 1.0, 1.0, 5.0, 1.0, 10.0, 6.0, 4.0, 10.0, 3.0, 6.0, 3.0, 3.0, 6.0, 8.0, 8.0, 5.0, 9.0, 3.0, 6.0, 5.0, 9.0, 8.0, 5.0, 6.0, 2.0, 6.0, 3.0, 5.0, 9.0, 4.0, 2.0, 10.0, 3.0, 10.0, 5.0, 2.0, 5.0, 8.0, 3.0, 10.0, 3.0, 2.0, 1.0, 8.0, 4.0, 8.0, 8.0, 7.0, 9.0, 9.0]
global b_y = 10
global p = [0.334, 0.273, 0.363, 0.901, 0.001, 0.946, 0.956, 0.47, 0.784, 0.992, 0.77, 0.362, 0.176, 0.142, 0.526, 0.335, 0.934, 0.214, 0.554, 0.938, 0.503, 0.686, 0.757, 0.137, 0.001, 0.731, 0.177, 0.346, 0.176, 0.502, 0.965, 0.098, 0.955, 0.807, 0.017, 0.365, 0.682, 0.665, 0.073, 0.373, 0.333, 0.074, 0.939, 0.315, 0.559, 0.719, 0.701, 0.036, 0.773, 0.248, 0.828, 0.029, 0.91, 0.599, 0.642, 0.429, 0.274, 0.121, 0.866, 0.852, 0.675, 0.875, 0.477, 0.938, 0.958, 0.517, 0.663, 0.26, 0.236, 0.927, 0.405, 0.725, 0.603, 0.041, 0.295, 0.19, 0.567, 0.14, 0.177, 0.553, 0.24, 0.609, 0.669, 0.471, 0.328, 0.788, 0.105, 0.782, 0.235, 0.647, 0.762, 0.255, 0.351, 0.052, 0.364, 0.789, 0.992, 0.627, 0.354, 0.056, 0.593, 0.477, 0.576, 0.639, 0.024, 0.93, 0.011, 0.476, 0.113, 0.272, 0.666, 0.807, 0.518, 0.593, 0.81, 0.094, 0.455, 0.831, 0.296, 0.24, 0.366, 0.635, 0.983, 0.318, 0.372, 0.984, 0.621, 0.469, 0.85, 0.28, 0.834, 0.245, 0.49, 0.444, 0.865, 0.458, 0.056, 0.647, 0.517, 0.684, 0.776, 0.607, 0.77, 0.193, 0.9, 0.052, 0.073, 0.34, 0.481, 0.192, 0.768, 0.289, 0.71, 0.306, 0.826, 0.58, 0.245, 0.303, 0.951, 0.313, 0.235, 0.86, 0.882, 0.446, 0.967, 0.968, 0.12, 0.781, 0.17, 0.069, 0.377, 0.336, 0.664, 0.671, 0.542, 0.526, 0.146, 0.671, 0.388, 0.959, 0.77, 0.355, 0.437, 0.846, 0.02, 0.691, 0.075, 0.633, 0.686, 0.843, 0.909, 0.502, 0.178, 0.191, 0.66, 0.171, 0.91, 0.847, 0.299, 0.016, 0.158, 0.171, 0.302, 0.334, 0.249, 0.082, 0.284, 0.135, 0.613, 0.015, 0.059, 0.159, 0.929, 0.157, 0.946, 0.875, 0.361, 0.119, 0.821, 0.861, 0.322, 0.601, 0.861, 0.01, 0.807, 0.997, 0.212, 0.507, 0.305, 0.11, 0.878, 0.698, 0.872, 0.504, 0.005, 0.689, 0.872, 0.058]
global q = [0.421, 0.323, 0.789, 0.96, 0.559, 0.989, 0.964, 0.965, 0.882, 0.992, 0.953, 0.955, 0.622, 0.485, 0.673, 0.613, 0.991, 0.47, 0.616, 0.954, 0.713, 0.86, 0.77, 0.367, 0.296, 0.959, 0.315, 0.793, 0.452, 0.686, 0.982, 0.173, 0.999, 0.831, 0.348, 0.375, 0.709, 0.73, 0.549, 0.903, 0.341, 0.912, 0.999, 0.935, 0.836, 0.761, 0.785, 0.201, 0.915, 0.678, 0.849, 0.493, 0.969, 0.632, 0.744, 0.486, 0.328, 0.232, 0.969, 0.9, 0.853, 0.947, 0.752, 0.975, 0.975, 0.944, 0.762, 0.702, 0.763, 0.991, 0.916, 0.835, 0.96, 0.734, 0.58, 0.581, 0.742, 0.611, 0.561, 0.683, 0.561, 0.809, 0.757, 0.937, 0.588, 0.856, 0.653, 0.948, 0.697, 0.777, 0.993, 0.983, 0.768, 0.343, 0.663, 0.796, 0.993, 0.902, 0.867, 0.916, 0.793, 0.782, 0.825, 0.684, 0.485, 0.956, 0.86, 0.608, 0.61, 0.804, 0.925, 0.919, 0.557, 0.823, 0.99, 0.98, 0.758, 0.834, 0.95, 0.592, 0.663, 0.901, 0.988, 0.978, 0.995, 0.999, 0.988, 0.925, 0.908, 0.72, 0.959, 0.32, 0.707, 0.921, 0.949, 0.878, 0.744, 0.66, 0.84, 0.722, 0.892, 0.892, 0.802, 0.747, 0.984, 0.529, 0.318, 0.809, 0.645, 0.98, 0.794, 0.748, 0.769, 0.529, 0.874, 0.903, 0.33, 0.877, 0.995, 0.339, 0.424, 0.922, 0.949, 0.768, 0.991, 0.989, 0.606, 0.989, 0.665, 0.61, 0.54, 0.889, 0.826, 0.685, 0.774, 0.705, 0.42, 0.86, 0.51, 0.969, 0.919, 0.508, 0.471, 0.862, 0.077, 0.86, 0.882, 0.662, 0.93, 0.99, 0.964, 0.859, 0.564, 0.862, 0.984, 0.782, 0.916, 0.956, 0.593, 0.998, 0.854, 0.211, 0.433, 0.974, 0.729, 0.683, 0.648, 0.374, 0.718, 0.166, 0.172, 0.963, 0.987, 0.92, 0.978, 0.92, 0.463, 0.875, 0.885, 0.999, 0.689, 0.609, 0.968, 0.782, 0.947, 0.999, 0.257, 0.762, 0.455, 0.689, 0.892, 0.815, 0.956, 0.52, 0.759, 0.807, 0.946, 0.399]
global origin = 1
global destination = 50