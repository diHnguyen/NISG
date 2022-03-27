global arcs = [1 14; 1 22; 1 28; 1 38; 2 6; 2 10; 2 12; 2 31; 3 4; 3 22; 3 42; 3 44; 4 9; 4 10; 4 16; 4 27; 4 36; 4 44; 5 17; 5 21; 5 29; 5 30; 5 32; 5 35; 5 36; 5 50; 6 14; 6 17; 6 26; 6 45; 6 49; 7 16; 7 17; 7 26; 7 28; 7 35; 7 36; 7 41; 8 12; 8 44; 8 47; 8 48; 9 4; 9 6; 9 44; 10 6; 10 29; 10 44; 11 24; 11 41; 12 6; 12 30; 12 38; 12 44; 13 2; 13 3; 13 17; 13 21; 13 32; 13 42; 13 45; 13 46; 14 2; 14 5; 14 9; 14 18; 14 22; 14 23; 14 38; 15 12; 15 16; 15 18; 15 21; 15 23; 15 25; 15 30; 15 31; 15 33; 15 34; 15 48; 16 29; 16 33; 16 44; 16 45; 16 46; 17 14; 17 27; 17 45; 17 48; 18 14; 18 25; 18 26; 18 30; 18 33; 18 34; 18 35; 18 41; 18 45; 18 49; 19 7; 19 32; 20 27; 20 42; 21 17; 21 30; 21 47; 22 5; 22 9; 22 13; 22 16; 22 28; 22 46; 23 3; 23 8; 23 20; 23 29; 23 37; 24 7; 24 18; 24 20; 24 23; 24 36; 24 45; 25 3; 25 19; 25 23; 25 24; 26 4; 26 6; 26 18; 26 23; 26 30; 26 34; 26 37; 26 38; 26 39; 26 50; 27 7; 27 23; 27 33; 27 35; 27 49; 28 4; 28 6; 28 8; 28 13; 28 30; 28 34; 28 37; 28 42; 28 48; 28 50; 29 2; 29 19; 29 42; 29 48; 30 6; 30 14; 30 19; 30 28; 30 29; 30 37; 30 45; 30 46; 31 10; 31 11; 32 5; 32 19; 33 7; 33 15; 33 18; 33 29; 33 46; 33 47; 34 10; 34 12; 34 16; 34 28; 34 32; 34 50; 35 13; 36 14; 36 22; 36 23; 36 27; 36 30; 36 35; 36 39; 36 40; 36 47; 37 24; 37 28; 37 40; 37 47; 38 5; 38 21; 39 25; 39 43; 39 44; 40 10; 40 17; 40 30; 40 34; 40 42; 41 6; 41 7; 41 10; 41 14; 41 37; 41 38; 41 43; 42 11; 42 19; 42 24; 42 35; 42 48; 42 49; 42 50; 43 14; 43 17; 43 26; 43 28; 43 29; 44 19; 44 23; 44 31; 44 34; 44 37; 45 11; 45 13; 46 2; 46 10; 46 16; 46 18; 46 21; 46 27; 46 37; 47 2; 47 6; 47 14; 47 17; 47 29; 47 38; 48 9; 48 35; 49 24; 49 30; 49 33; 49 42]
global d_x = [2.0, 2.0, 9.0, 10.0, 7.0, 10.0, 10.0, 6.0, 3.0, 8.0, 5.0, 10.0, 10.0, 10.0, 7.0, 8.0, 1.0, 3.0, 8.0, 1.0, 6.0, 7.0, 2.0, 4.0, 5.0, 3.0, 5.0, 9.0, 2.0, 4.0, 5.0, 1.0, 5.0, 9.0, 10.0, 10.0, 3.0, 10.0, 3.0, 3.0, 10.0, 7.0, 7.0, 9.0, 2.0, 9.0, 4.0, 8.0, 10.0, 4.0, 6.0, 8.0, 9.0, 7.0, 2.0, 5.0, 5.0, 9.0, 9.0, 8.0, 10.0, 10.0, 6.0, 1.0, 5.0, 7.0, 1.0, 1.0, 5.0, 3.0, 8.0, 2.0, 6.0, 4.0, 6.0, 6.0, 2.0, 7.0, 10.0, 9.0, 7.0, 2.0, 3.0, 6.0, 7.0, 5.0, 1.0, 7.0, 4.0, 8.0, 2.0, 3.0, 9.0, 9.0, 4.0, 2.0, 9.0, 2.0, 10.0, 1.0, 5.0, 3.0, 9.0, 9.0, 2.0, 10.0, 9.0, 1.0, 8.0, 6.0, 2.0, 4.0, 3.0, 5.0, 5.0, 7.0, 9.0, 3.0, 1.0, 7.0, 1.0, 9.0, 8.0, 5.0, 2.0, 2.0, 1.0, 6.0, 4.0, 4.0, 7.0, 7.0, 3.0, 1.0, 10.0, 2.0, 8.0, 2.0, 4.0, 7.0, 7.0, 8.0, 3.0, 10.0, 7.0, 9.0, 3.0, 1.0, 5.0, 2.0, 8.0, 9.0, 7.0, 3.0, 8.0, 9.0, 10.0, 7.0, 8.0, 4.0, 4.0, 5.0, 6.0, 1.0, 7.0, 5.0, 4.0, 1.0, 6.0, 8.0, 6.0, 4.0, 4.0, 1.0, 8.0, 1.0, 6.0, 10.0, 10.0, 9.0, 7.0, 3.0, 8.0, 10.0, 1.0, 8.0, 7.0, 4.0, 5.0, 3.0, 1.0, 10.0, 1.0, 1.0, 1.0, 3.0, 10.0, 5.0, 6.0, 2.0, 3.0, 9.0, 3.0, 2.0, 1.0, 4.0, 4.0, 3.0, 4.0, 9.0, 6.0, 9.0, 7.0, 6.0, 4.0, 7.0, 1.0, 6.0, 6.0, 3.0, 8.0, 6.0, 4.0, 1.0, 6.0, 7.0, 7.0, 2.0, 7.0, 7.0, 8.0, 4.0, 6.0, 2.0, 5.0, 8.0, 8.0, 3.0, 5.0, 2.0, 10.0, 2.0, 4.0, 2.0, 5.0, 2.0, 8.0, 1.0, 8.0]
global b_x = 5
global d_y = [9.0, 3.0, 2.0, 5.0, 2.0, 5.0, 1.0, 2.0, 5.0, 10.0, 1.0, 3.0, 9.0, 9.0, 7.0, 4.0, 6.0, 3.0, 3.0, 8.0, 2.0, 1.0, 10.0, 4.0, 10.0, 7.0, 4.0, 3.0, 10.0, 3.0, 10.0, 2.0, 5.0, 7.0, 6.0, 4.0, 10.0, 2.0, 6.0, 10.0, 9.0, 3.0, 10.0, 6.0, 10.0, 3.0, 3.0, 3.0, 4.0, 5.0, 1.0, 8.0, 2.0, 9.0, 2.0, 9.0, 3.0, 1.0, 9.0, 9.0, 4.0, 8.0, 8.0, 1.0, 7.0, 10.0, 8.0, 3.0, 1.0, 9.0, 4.0, 2.0, 2.0, 1.0, 10.0, 9.0, 6.0, 6.0, 5.0, 9.0, 5.0, 2.0, 1.0, 4.0, 6.0, 10.0, 4.0, 8.0, 7.0, 8.0, 3.0, 5.0, 2.0, 8.0, 3.0, 7.0, 1.0, 6.0, 5.0, 3.0, 2.0, 1.0, 3.0, 8.0, 8.0, 5.0, 6.0, 6.0, 7.0, 10.0, 2.0, 4.0, 2.0, 1.0, 5.0, 3.0, 7.0, 1.0, 3.0, 7.0, 5.0, 4.0, 1.0, 4.0, 8.0, 2.0, 2.0, 9.0, 6.0, 6.0, 3.0, 9.0, 9.0, 3.0, 8.0, 8.0, 4.0, 3.0, 3.0, 7.0, 10.0, 1.0, 9.0, 1.0, 10.0, 10.0, 4.0, 8.0, 4.0, 7.0, 1.0, 7.0, 5.0, 3.0, 3.0, 7.0, 9.0, 3.0, 8.0, 3.0, 7.0, 4.0, 1.0, 2.0, 10.0, 5.0, 8.0, 9.0, 4.0, 6.0, 5.0, 4.0, 4.0, 9.0, 3.0, 6.0, 10.0, 6.0, 7.0, 3.0, 10.0, 8.0, 8.0, 3.0, 1.0, 9.0, 9.0, 7.0, 4.0, 8.0, 8.0, 6.0, 1.0, 8.0, 8.0, 3.0, 7.0, 3.0, 1.0, 9.0, 3.0, 2.0, 3.0, 4.0, 1.0, 9.0, 1.0, 7.0, 3.0, 7.0, 10.0, 1.0, 5.0, 4.0, 10.0, 10.0, 4.0, 6.0, 7.0, 9.0, 3.0, 1.0, 1.0, 4.0, 2.0, 8.0, 8.0, 2.0, 6.0, 5.0, 3.0, 5.0, 5.0, 8.0, 4.0, 5.0, 3.0, 10.0, 4.0, 10.0, 3.0, 10.0, 5.0, 7.0, 6.0, 1.0, 3.0, 10.0, 7.0]
global b_y = 10
global p = [0.4, 0.197, 0.114, 0.904, 0.888, 0.159, 0.896, 0.98, 0.56, 0.627, 0.088, 0.595, 0.787, 0.3, 0.536, 0.496, 0.616, 0.816, 0.118, 0.138, 0.207, 0.145, 0.597, 0.969, 0.122, 0.835, 0.492, 0.832, 0.627, 0.099, 0.949, 0.702, 0.611, 0.998, 0.904, 0.549, 0.248, 0.555, 0.788, 0.083, 0.899, 0.895, 0.504, 0.764, 0.884, 0.145, 0.033, 0.522, 0.33, 0.593, 0.274, 0.957, 0.032, 0.166, 0.112, 0.744, 0.648, 0.536, 0.407, 0.827, 0.681, 0.519, 0.284, 0.046, 0.426, 0.047, 0.286, 0.337, 0.632, 0.878, 0.834, 0.134, 0.661, 0.198, 0.767, 0.455, 0.132, 0.302, 0.512, 0.621, 0.085, 0.39, 0.135, 0.379, 0.407, 0.962, 0.506, 0.522, 0.823, 0.828, 0.671, 0.132, 0.075, 0.369, 0.975, 0.878, 0.562, 0.769, 0.037, 0.886, 0.192, 0.057, 0.898, 0.547, 0.181, 0.779, 0.929, 0.336, 0.616, 0.209, 0.265, 0.066, 0.096, 0.971, 0.013, 0.832, 0.756, 0.824, 0.974, 0.481, 0.749, 0.715, 0.58, 0.021, 0.229, 0.652, 0.03, 0.074, 0.881, 0.206, 0.162, 0.776, 0.86, 0.425, 0.58, 0.307, 0.916, 0.646, 0.359, 0.181, 0.689, 0.34, 0.72, 0.987, 0.807, 0.628, 0.789, 0.539, 0.961, 0.246, 0.352, 0.653, 0.931, 0.346, 0.984, 0.97, 0.933, 0.149, 0.395, 0.099, 0.586, 0.515, 0.95, 0.444, 0.746, 0.109, 0.757, 0.804, 0.321, 0.63, 0.718, 0.284, 0.805, 0.372, 0.38, 0.455, 0.955, 0.07, 0.216, 0.676, 0.7, 0.88, 0.279, 0.898, 0.282, 0.917, 0.216, 0.436, 0.434, 0.345, 0.772, 0.533, 0.029, 0.827, 0.472, 0.508, 0.532, 0.376, 0.539, 0.468, 0.822, 0.383, 0.931, 0.118, 0.649, 0.609, 0.336, 0.548, 0.178, 0.344, 0.44, 0.743, 0.775, 0.394, 0.588, 0.342, 0.774, 0.209, 0.353, 0.919, 0.291, 0.217, 0.456, 0.437, 0.496, 0.961, 0.273, 0.429, 0.039, 0.021, 0.563, 0.486, 0.302, 0.346, 0.455, 0.656, 0.02, 0.435, 0.71, 0.209, 0.549, 0.919, 0.891, 0.187, 0.575, 0.066, 0.296, 0.675, 0.225]
global q = [0.799, 0.703, 0.149, 0.962, 0.921, 0.489, 0.996, 0.99, 0.969, 0.876, 0.123, 0.637, 0.933, 0.417, 0.732, 0.788, 0.883, 0.835, 0.481, 0.149, 0.409, 0.522, 0.652, 0.978, 0.461, 0.839, 0.798, 0.997, 0.802, 0.175, 0.951, 0.842, 0.74, 0.998, 0.934, 0.975, 0.656, 0.978, 0.982, 0.59, 0.958, 0.989, 0.935, 0.999, 0.968, 0.258, 0.489, 0.898, 0.753, 0.913, 0.998, 0.97, 0.357, 0.87, 0.349, 0.991, 0.853, 0.593, 0.659, 0.873, 0.996, 0.966, 0.287, 0.614, 0.614, 0.817, 0.938, 0.776, 0.636, 0.911, 0.839, 0.763, 0.683, 0.358, 0.877, 0.69, 0.283, 0.309, 0.733, 0.877, 0.276, 0.616, 0.705, 0.856, 0.81, 0.969, 0.777, 0.903, 0.938, 0.859, 0.868, 0.998, 0.489, 0.868, 0.999, 0.903, 0.896, 0.807, 0.665, 0.968, 0.502, 0.1, 0.961, 0.788, 0.348, 0.817, 0.93, 0.452, 0.845, 0.689, 0.31, 0.558, 0.955, 0.996, 0.282, 0.91, 0.932, 0.854, 0.984, 0.764, 0.849, 0.811, 0.818, 0.879, 0.625, 0.826, 0.799, 0.295, 0.922, 0.57, 0.553, 0.945, 0.982, 0.471, 0.925, 0.59, 0.938, 0.967, 0.519, 0.889, 0.854, 0.679, 0.813, 0.993, 0.918, 0.661, 0.907, 0.728, 0.969, 0.932, 0.453, 0.953, 0.959, 0.991, 0.997, 0.977, 0.961, 0.522, 0.739, 0.734, 0.685, 0.869, 0.971, 0.514, 0.845, 0.851, 0.851, 0.965, 0.988, 0.815, 0.815, 0.923, 0.987, 0.665, 0.577, 0.646, 0.971, 0.622, 0.668, 0.855, 0.991, 0.881, 0.954, 0.954, 0.666, 0.955, 0.698, 0.982, 0.981, 0.533, 0.807, 0.557, 0.783, 0.835, 0.534, 0.791, 0.627, 0.877, 0.634, 0.491, 0.89, 0.511, 0.956, 0.258, 0.82, 0.721, 0.637, 0.996, 0.679, 0.585, 0.768, 0.896, 0.838, 0.89, 0.593, 0.35, 0.944, 0.964, 0.625, 0.995, 0.863, 0.931, 0.458, 0.64, 0.71, 0.992, 0.771, 0.685, 0.27, 0.956, 0.897, 0.963, 0.806, 0.931, 0.739, 0.858, 0.475, 0.801, 0.73, 0.722, 0.623, 0.94, 0.901, 0.31, 0.671, 0.424, 0.721, 0.936, 0.38]
global origin = 1
global destination = 50