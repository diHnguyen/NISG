global arcs = [1 18; 1 44; 1 49; 2 20; 2 37; 2 47; 3 20; 3 37; 3 39; 3 40; 3 41; 3 42; 4 3; 4 13; 4 39; 5 27; 5 33; 6 15; 6 22; 6 31; 6 37; 6 40; 6 42; 6 48; 6 50; 7 12; 7 20; 7 28; 7 31; 7 35; 7 50; 8 16; 8 30; 9 3; 9 15; 9 21; 9 23; 9 27; 9 36; 10 19; 10 32; 10 37; 10 43; 11 8; 11 34; 11 39; 12 3; 12 14; 12 41; 12 45; 13 32; 13 39; 13 45; 14 3; 14 6; 14 33; 14 41; 14 46; 14 47; 14 48; 15 2; 15 3; 15 12; 15 31; 15 34; 15 36; 16 6; 16 27; 16 29; 16 32; 16 47; 16 50; 17 5; 17 12; 17 20; 17 28; 17 36; 17 39; 17 46; 18 7; 18 8; 18 28; 18 49; 19 23; 19 34; 20 7; 20 8; 20 12; 20 25; 20 26; 20 30; 20 37; 20 40; 20 43; 21 9; 21 13; 21 20; 21 23; 22 3; 22 6; 22 13; 22 23; 22 25; 22 37; 22 41; 22 45; 22 46; 22 48; 23 5; 23 15; 23 19; 23 35; 23 49; 24 18; 24 33; 24 34; 24 44; 25 2; 25 11; 25 26; 25 32; 26 27; 26 33; 26 41; 26 45; 27 3; 27 8; 27 50; 28 5; 28 29; 28 41; 28 42; 29 4; 29 8; 29 9; 29 18; 29 25; 29 40; 29 42; 29 45; 30 13; 30 20; 30 24; 30 42; 30 47; 31 3; 31 20; 31 25; 32 13; 32 27; 32 35; 32 38; 32 39; 33 15; 33 37; 34 3; 34 27; 34 42; 34 44; 35 4; 35 13; 35 15; 35 17; 35 26; 35 29; 35 42; 36 15; 36 44; 36 45; 36 49; 36 50; 37 13; 37 26; 37 41; 37 44; 38 3; 38 9; 38 14; 38 19; 38 25; 38 29; 39 14; 39 19; 39 26; 39 35; 39 40; 39 46; 39 48; 40 3; 40 11; 40 14; 40 41; 41 8; 41 11; 41 13; 41 18; 41 19; 41 24; 41 29; 41 31; 42 2; 42 10; 42 26; 42 31; 42 33; 42 34; 43 12; 43 19; 43 25; 43 30; 43 34; 43 46; 44 20; 44 38; 44 43; 44 47; 44 48; 45 4; 45 42; 46 2; 46 10; 46 21; 46 32; 46 34; 47 6; 47 7; 47 27; 47 44; 48 4; 48 8; 48 13; 48 15; 48 33; 48 42; 48 45; 48 50; 49 23; 49 28]
global d_x = [2.0, 7.0, 3.0, 4.0, 3.0, 8.0, 5.0, 1.0, 8.0, 2.0, 7.0, 5.0, 7.0, 10.0, 6.0, 10.0, 10.0, 4.0, 5.0, 7.0, 8.0, 3.0, 1.0, 4.0, 6.0, 1.0, 3.0, 4.0, 4.0, 5.0, 10.0, 1.0, 10.0, 6.0, 5.0, 5.0, 5.0, 9.0, 3.0, 1.0, 2.0, 6.0, 7.0, 7.0, 4.0, 4.0, 4.0, 7.0, 8.0, 4.0, 1.0, 2.0, 9.0, 3.0, 2.0, 4.0, 6.0, 1.0, 7.0, 2.0, 10.0, 4.0, 8.0, 5.0, 1.0, 2.0, 5.0, 9.0, 7.0, 2.0, 4.0, 1.0, 8.0, 1.0, 1.0, 6.0, 9.0, 2.0, 9.0, 3.0, 9.0, 9.0, 1.0, 9.0, 7.0, 2.0, 3.0, 7.0, 1.0, 3.0, 1.0, 4.0, 7.0, 9.0, 4.0, 9.0, 9.0, 6.0, 8.0, 7.0, 1.0, 6.0, 5.0, 3.0, 2.0, 5.0, 10.0, 7.0, 8.0, 4.0, 9.0, 6.0, 3.0, 3.0, 4.0, 1.0, 3.0, 2.0, 8.0, 1.0, 9.0, 6.0, 2.0, 4.0, 9.0, 4.0, 2.0, 3.0, 5.0, 8.0, 3.0, 5.0, 7.0, 10.0, 10.0, 9.0, 8.0, 7.0, 3.0, 10.0, 3.0, 5.0, 3.0, 1.0, 4.0, 1.0, 7.0, 5.0, 6.0, 1.0, 7.0, 6.0, 8.0, 1.0, 3.0, 2.0, 3.0, 5.0, 9.0, 10.0, 6.0, 6.0, 2.0, 2.0, 2.0, 8.0, 1.0, 7.0, 5.0, 9.0, 8.0, 3.0, 2.0, 8.0, 6.0, 7.0, 8.0, 2.0, 2.0, 8.0, 9.0, 10.0, 7.0, 9.0, 8.0, 5.0, 1.0, 2.0, 9.0, 5.0, 2.0, 9.0, 7.0, 3.0, 2.0, 6.0, 1.0, 10.0, 6.0, 8.0, 9.0, 3.0, 4.0, 9.0, 4.0, 7.0, 4.0, 10.0, 5.0, 5.0, 1.0, 2.0, 8.0, 8.0, 5.0, 9.0, 6.0, 1.0, 1.0, 3.0, 9.0, 2.0, 3.0, 6.0, 7.0, 9.0, 8.0, 10.0, 10.0, 6.0, 6.0, 9.0, 10.0, 8.0, 8.0, 10.0, 1.0, 1.0]
global b_x = 5
global d_y = [10.0, 1.0, 2.0, 1.0, 5.0, 1.0, 10.0, 3.0, 6.0, 5.0, 8.0, 7.0, 6.0, 9.0, 8.0, 7.0, 6.0, 6.0, 6.0, 9.0, 1.0, 4.0, 1.0, 1.0, 9.0, 5.0, 3.0, 4.0, 2.0, 7.0, 7.0, 2.0, 5.0, 4.0, 4.0, 2.0, 10.0, 5.0, 6.0, 4.0, 1.0, 2.0, 5.0, 5.0, 9.0, 1.0, 7.0, 4.0, 2.0, 9.0, 9.0, 10.0, 1.0, 4.0, 5.0, 6.0, 9.0, 7.0, 1.0, 9.0, 1.0, 2.0, 1.0, 10.0, 7.0, 3.0, 10.0, 5.0, 10.0, 8.0, 6.0, 10.0, 5.0, 7.0, 6.0, 4.0, 5.0, 8.0, 5.0, 8.0, 4.0, 2.0, 7.0, 9.0, 5.0, 9.0, 3.0, 2.0, 8.0, 5.0, 10.0, 1.0, 10.0, 7.0, 2.0, 5.0, 6.0, 1.0, 5.0, 4.0, 8.0, 9.0, 3.0, 3.0, 3.0, 8.0, 9.0, 9.0, 1.0, 5.0, 1.0, 8.0, 4.0, 6.0, 9.0, 4.0, 5.0, 6.0, 2.0, 4.0, 6.0, 3.0, 7.0, 10.0, 7.0, 4.0, 3.0, 4.0, 1.0, 2.0, 8.0, 8.0, 7.0, 5.0, 10.0, 9.0, 5.0, 6.0, 2.0, 6.0, 10.0, 7.0, 1.0, 10.0, 3.0, 6.0, 2.0, 10.0, 2.0, 7.0, 2.0, 7.0, 1.0, 8.0, 1.0, 6.0, 4.0, 8.0, 1.0, 3.0, 7.0, 3.0, 1.0, 8.0, 4.0, 4.0, 10.0, 2.0, 2.0, 7.0, 4.0, 10.0, 4.0, 3.0, 10.0, 8.0, 4.0, 10.0, 10.0, 7.0, 2.0, 8.0, 4.0, 3.0, 1.0, 10.0, 8.0, 4.0, 1.0, 5.0, 8.0, 5.0, 6.0, 2.0, 2.0, 10.0, 10.0, 6.0, 5.0, 6.0, 7.0, 4.0, 10.0, 6.0, 10.0, 7.0, 9.0, 4.0, 3.0, 6.0, 10.0, 2.0, 5.0, 7.0, 6.0, 3.0, 7.0, 4.0, 10.0, 10.0, 6.0, 2.0, 5.0, 8.0, 6.0, 7.0, 6.0, 5.0, 8.0, 4.0, 8.0, 1.0, 4.0, 3.0, 3.0, 10.0, 4.0, 10.0]
global b_y = 10
global p = [0.709, 0.051, 0.49, 0.896, 0.368, 0.294, 0.598, 0.897, 0.677, 0.786, 0.601, 0.205, 0.578, 0.785, 0.958, 0.51, 0.689, 0.131, 0.943, 0.766, 0.998, 0.6, 0.7, 0.133, 0.239, 0.963, 0.998, 0.165, 0.529, 0.756, 0.651, 0.777, 0.792, 0.966, 0.202, 0.028, 0.629, 0.352, 0.285, 0.921, 0.665, 0.818, 0.981, 0.301, 0.514, 0.883, 0.934, 0.817, 0.725, 0.993, 0.972, 0.327, 0.555, 0.575, 0.478, 0.042, 0.356, 0.381, 0.025, 0.018, 0.659, 0.727, 0.223, 0.279, 0.802, 0.461, 0.503, 0.902, 0.605, 0.479, 0.295, 0.011, 0.348, 0.007, 0.24, 0.973, 0.241, 0.807, 0.24, 0.472, 0.41, 0.467, 0.558, 0.644, 0.672, 0.105, 0.215, 0.939, 0.113, 0.858, 0.727, 0.507, 0.77, 0.657, 0.487, 0.784, 0.26, 0.603, 0.781, 0.803, 0.374, 0.274, 0.643, 0.307, 0.44, 0.386, 0.286, 0.632, 0.521, 0.161, 0.789, 0.51, 0.427, 0.502, 0.138, 0.036, 0.186, 0.016, 0.594, 0.883, 0.806, 0.143, 0.217, 0.265, 0.689, 0.077, 0.835, 0.491, 0.951, 0.658, 0.304, 0.749, 0.828, 0.782, 0.581, 0.229, 0.98, 0.131, 0.038, 0.342, 0.67, 0.443, 0.532, 0.795, 0.148, 0.286, 0.824, 0.649, 0.051, 0.464, 0.599, 0.435, 0.054, 0.829, 0.461, 0.257, 0.545, 0.104, 0.534, 0.971, 0.818, 0.463, 0.303, 0.774, 0.774, 0.8, 0.247, 0.971, 0.184, 0.743, 0.59, 0.871, 0.621, 0.132, 0.815, 0.448, 0.711, 0.896, 0.986, 0.525, 0.295, 0.798, 0.596, 0.143, 0.445, 0.723, 0.286, 0.992, 0.57, 0.845, 0.822, 0.035, 0.567, 0.683, 0.035, 0.359, 0.248, 0.64, 0.659, 0.499, 0.403, 0.112, 0.103, 0.943, 0.642, 0.668, 0.774, 0.734, 0.786, 0.606, 0.679, 0.436, 0.968, 0.971, 0.046, 0.099, 0.56, 0.897, 0.404, 0.116, 0.452, 0.133, 0.9, 0.329, 0.401, 0.734, 0.525, 0.809, 0.516, 0.339, 0.072, 0.015, 0.848, 0.012, 0.169, 0.718, 0.227, 0.962]
global q = [0.936, 0.327, 0.999, 0.953, 0.93, 0.971, 0.814, 0.923, 0.691, 0.96, 0.876, 0.467, 0.958, 0.904, 0.962, 0.58, 0.835, 0.499, 0.951, 0.769, 0.999, 0.825, 0.762, 0.645, 0.665, 0.987, 0.998, 0.92, 0.726, 0.836, 0.858, 0.863, 0.917, 0.977, 0.229, 0.276, 0.634, 0.668, 0.841, 0.922, 0.696, 0.861, 0.995, 0.64, 0.73, 0.963, 0.945, 0.988, 0.859, 0.995, 0.984, 0.929, 0.722, 0.775, 0.67, 0.138, 0.457, 0.81, 0.028, 0.725, 0.868, 0.802, 0.999, 0.445, 0.85, 0.552, 0.757, 0.924, 0.737, 0.701, 0.626, 0.2, 0.937, 0.138, 0.828, 0.985, 0.277, 0.96, 0.628, 0.725, 0.469, 0.538, 0.682, 0.757, 0.69, 0.402, 0.387, 0.967, 0.464, 0.89, 0.901, 0.866, 0.987, 0.982, 0.534, 0.887, 0.775, 0.947, 0.794, 0.828, 0.971, 0.407, 0.986, 0.597, 0.985, 0.513, 0.905, 0.936, 0.709, 0.535, 0.863, 0.594, 0.61, 0.854, 0.851, 0.487, 0.604, 0.434, 0.775, 0.951, 0.89, 0.903, 0.238, 0.385, 0.913, 0.654, 0.958, 0.932, 0.968, 0.681, 0.452, 0.9, 0.928, 0.878, 0.962, 0.536, 0.995, 0.89, 0.278, 0.66, 0.769, 0.476, 0.601, 0.817, 0.915, 0.561, 0.849, 0.7, 0.605, 0.497, 0.661, 0.607, 0.529, 0.925, 0.715, 0.534, 0.861, 0.356, 0.918, 0.973, 0.892, 0.677, 0.398, 0.841, 0.852, 0.923, 0.848, 0.982, 0.552, 0.91, 0.788, 0.997, 0.906, 0.708, 0.849, 0.731, 0.802, 0.942, 0.989, 0.906, 0.406, 0.917, 0.997, 0.936, 0.567, 0.992, 0.959, 0.997, 0.691, 0.929, 0.912, 0.915, 0.614, 0.96, 0.27, 0.999, 0.748, 0.772, 0.946, 0.858, 0.963, 0.468, 0.265, 0.949, 0.731, 0.868, 0.949, 0.99, 0.972, 0.677, 0.902, 0.519, 0.989, 0.973, 0.674, 0.87, 0.655, 0.963, 0.517, 0.931, 0.986, 0.631, 0.992, 0.759, 0.548, 0.797, 0.682, 0.881, 0.577, 0.734, 0.515, 0.848, 0.876, 0.26, 0.207, 0.787, 0.422, 0.962]
global origin = 1
global destination = 50