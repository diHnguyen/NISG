global arcs = [1 10; 1 15; 1 22; 2 19; 2 26; 2 35; 2 40; 2 41; 3 4; 3 8; 3 28; 3 39; 3 40; 3 44; 4 6; 4 9; 4 10; 4 12; 4 25; 4 37; 4 41; 4 43; 5 20; 5 23; 6 13; 6 15; 6 17; 6 20; 6 22; 6 40; 7 4; 7 13; 7 21; 7 24; 7 34; 7 39; 7 46; 8 9; 8 10; 8 20; 8 27; 8 34; 8 40; 8 44; 8 46; 8 47; 9 4; 9 12; 9 35; 9 40; 9 50; 10 6; 10 9; 10 39; 10 43; 10 48; 11 3; 11 17; 11 46; 11 47; 12 13; 12 18; 12 21; 12 38; 13 15; 13 19; 13 36; 13 41; 13 44; 13 50; 14 37; 14 39; 14 44; 15 21; 15 34; 16 11; 16 12; 16 19; 16 32; 16 33; 16 45; 16 50; 17 10; 17 18; 17 21; 17 45; 18 4; 18 14; 18 15; 18 29; 18 33; 18 46; 19 15; 19 22; 19 24; 19 45; 19 49; 20 7; 20 8; 20 19; 20 30; 20 33; 20 46; 21 3; 21 18; 21 19; 21 28; 21 32; 22 6; 22 13; 22 16; 22 29; 22 50; 23 5; 23 35; 23 36; 23 47; 24 20; 24 29; 25 11; 25 16; 25 36; 25 40; 25 43; 26 15; 26 27; 27 2; 27 24; 27 26; 27 34; 27 40; 28 11; 28 44; 28 45; 29 2; 29 42; 29 43; 29 44; 29 48; 30 7; 30 21; 30 24; 30 26; 30 27; 30 34; 31 2; 31 17; 31 18; 31 24; 31 27; 31 37; 32 11; 32 48; 33 6; 33 18; 33 24; 33 32; 33 43; 33 48; 34 20; 34 21; 34 26; 34 28; 34 30; 34 46; 35 4; 35 17; 35 20; 35 22; 35 32; 35 38; 35 50; 36 21; 36 23; 36 31; 36 49; 37 10; 37 32; 37 39; 37 42; 38 7; 38 25; 38 42; 38 44; 38 45; 38 46; 39 6; 39 21; 39 30; 40 7; 40 20; 40 27; 40 37; 40 45; 41 8; 41 17; 41 42; 41 46; 42 10; 42 15; 42 18; 42 19; 42 29; 42 32; 43 15; 43 32; 43 35; 43 50; 44 6; 44 15; 45 9; 45 15; 45 32; 46 23; 46 26; 46 33; 46 35; 46 38; 46 47; 47 23; 47 38; 47 48; 48 2; 48 6; 48 8; 48 24; 48 40; 49 2; 49 13; 49 16; 49 18; 49 24; 49 31]
global d_x = [5.0, 9.0, 8.0, 2.0, 6.0, 3.0, 1.0, 4.0, 1.0, 4.0, 7.0, 8.0, 1.0, 9.0, 4.0, 6.0, 4.0, 5.0, 6.0, 2.0, 4.0, 3.0, 1.0, 7.0, 8.0, 5.0, 3.0, 6.0, 10.0, 9.0, 6.0, 6.0, 10.0, 6.0, 10.0, 6.0, 7.0, 2.0, 4.0, 3.0, 2.0, 8.0, 6.0, 9.0, 2.0, 2.0, 4.0, 8.0, 3.0, 9.0, 10.0, 4.0, 10.0, 4.0, 5.0, 1.0, 3.0, 4.0, 7.0, 5.0, 4.0, 5.0, 3.0, 9.0, 6.0, 6.0, 7.0, 3.0, 10.0, 1.0, 1.0, 4.0, 2.0, 6.0, 9.0, 2.0, 5.0, 3.0, 10.0, 6.0, 4.0, 1.0, 9.0, 9.0, 8.0, 2.0, 5.0, 2.0, 4.0, 3.0, 3.0, 6.0, 10.0, 6.0, 8.0, 6.0, 4.0, 5.0, 6.0, 8.0, 9.0, 3.0, 9.0, 7.0, 3.0, 6.0, 6.0, 1.0, 5.0, 9.0, 9.0, 1.0, 10.0, 9.0, 7.0, 2.0, 3.0, 10.0, 2.0, 3.0, 7.0, 6.0, 6.0, 7.0, 3.0, 6.0, 10.0, 1.0, 1.0, 1.0, 9.0, 1.0, 1.0, 10.0, 8.0, 5.0, 9.0, 5.0, 6.0, 5.0, 5.0, 1.0, 1.0, 9.0, 10.0, 10.0, 10.0, 4.0, 10.0, 2.0, 9.0, 10.0, 4.0, 3.0, 4.0, 3.0, 5.0, 1.0, 5.0, 4.0, 4.0, 4.0, 1.0, 8.0, 10.0, 10.0, 3.0, 2.0, 6.0, 5.0, 10.0, 8.0, 6.0, 8.0, 10.0, 8.0, 4.0, 4.0, 10.0, 5.0, 7.0, 6.0, 2.0, 4.0, 9.0, 7.0, 9.0, 4.0, 9.0, 7.0, 4.0, 3.0, 4.0, 4.0, 4.0, 9.0, 3.0, 6.0, 7.0, 1.0, 1.0, 7.0, 3.0, 1.0, 10.0, 1.0, 10.0, 3.0, 7.0, 8.0, 3.0, 2.0, 3.0, 9.0, 9.0, 10.0, 4.0, 7.0, 5.0, 10.0, 9.0, 8.0, 2.0, 7.0, 9.0, 3.0, 7.0, 2.0, 10.0, 8.0, 6.0, 4.0, 8.0]
global b_x = 5
global d_y = [9.0, 3.0, 9.0, 10.0, 1.0, 4.0, 4.0, 2.0, 2.0, 6.0, 4.0, 3.0, 9.0, 8.0, 1.0, 5.0, 7.0, 4.0, 4.0, 7.0, 6.0, 3.0, 3.0, 6.0, 6.0, 10.0, 8.0, 9.0, 9.0, 6.0, 1.0, 4.0, 3.0, 3.0, 10.0, 4.0, 4.0, 9.0, 9.0, 1.0, 8.0, 2.0, 1.0, 1.0, 3.0, 7.0, 8.0, 9.0, 9.0, 1.0, 7.0, 2.0, 6.0, 2.0, 7.0, 5.0, 6.0, 7.0, 3.0, 10.0, 7.0, 7.0, 6.0, 8.0, 4.0, 7.0, 4.0, 2.0, 9.0, 9.0, 10.0, 1.0, 10.0, 2.0, 3.0, 10.0, 4.0, 6.0, 9.0, 3.0, 2.0, 1.0, 5.0, 6.0, 2.0, 7.0, 4.0, 5.0, 5.0, 2.0, 9.0, 8.0, 10.0, 1.0, 4.0, 7.0, 10.0, 6.0, 10.0, 5.0, 2.0, 6.0, 10.0, 4.0, 9.0, 4.0, 8.0, 5.0, 4.0, 6.0, 5.0, 8.0, 2.0, 6.0, 7.0, 5.0, 7.0, 8.0, 3.0, 1.0, 10.0, 5.0, 5.0, 5.0, 10.0, 8.0, 3.0, 3.0, 3.0, 7.0, 2.0, 9.0, 5.0, 8.0, 6.0, 7.0, 3.0, 6.0, 6.0, 4.0, 5.0, 3.0, 7.0, 9.0, 8.0, 7.0, 6.0, 8.0, 5.0, 9.0, 9.0, 2.0, 8.0, 8.0, 4.0, 8.0, 8.0, 7.0, 9.0, 3.0, 2.0, 8.0, 3.0, 8.0, 7.0, 10.0, 7.0, 7.0, 4.0, 4.0, 6.0, 1.0, 9.0, 8.0, 10.0, 2.0, 7.0, 2.0, 7.0, 1.0, 8.0, 5.0, 8.0, 2.0, 4.0, 4.0, 4.0, 7.0, 8.0, 6.0, 7.0, 8.0, 8.0, 9.0, 7.0, 8.0, 6.0, 4.0, 10.0, 9.0, 2.0, 9.0, 9.0, 3.0, 3.0, 2.0, 2.0, 7.0, 10.0, 1.0, 2.0, 5.0, 3.0, 8.0, 8.0, 6.0, 1.0, 4.0, 2.0, 5.0, 10.0, 10.0, 2.0, 5.0, 7.0, 3.0, 2.0, 8.0, 7.0, 3.0, 7.0, 1.0, 5.0]
global b_y = 10
global p = [0.472, 0.366, 0.587, 0.368, 0.239, 0.296, 0.259, 0.541, 0.114, 0.103, 0.404, 0.723, 0.44, 0.497, 0.487, 0.392, 0.829, 0.292, 0.583, 0.325, 0.702, 0.54, 0.27, 0.711, 0.298, 0.144, 0.786, 0.445, 0.541, 0.756, 0.763, 0.102, 0.861, 0.082, 0.374, 0.825, 0.368, 0.899, 0.864, 0.711, 0.375, 0.766, 0.923, 0.142, 0.977, 0.929, 0.636, 0.019, 0.718, 0.257, 0.383, 0.533, 0.959, 0.511, 0.514, 0.185, 0.4, 0.293, 0.961, 0.228, 0.658, 0.059, 0.989, 0.249, 0.517, 0.063, 0.351, 0.004, 0.74, 0.576, 0.283, 0.419, 0.295, 0.674, 0.97, 0.752, 0.681, 0.009, 0.007, 0.032, 0.099, 0.146, 0.608, 0.416, 0.548, 0.179, 0.717, 0.661, 0.559, 0.57, 0.687, 0.131, 0.275, 0.989, 0.705, 0.214, 0.254, 0.177, 0.644, 0.076, 0.924, 0.032, 0.97, 0.544, 0.869, 0.579, 0.544, 0.85, 0.337, 0.845, 0.932, 0.041, 0.337, 0.792, 0.32, 0.66, 0.946, 0.898, 0.595, 0.039, 0.008, 0.175, 0.723, 0.486, 0.848, 0.269, 0.49, 0.207, 0.368, 0.679, 0.882, 0.478, 0.837, 0.201, 0.491, 0.12, 0.019, 0.295, 0.099, 0.908, 0.601, 0.369, 0.914, 0.882, 0.481, 0.872, 0.797, 0.359, 0.268, 0.582, 0.047, 0.233, 0.974, 0.779, 0.298, 0.144, 0.761, 0.861, 0.182, 0.185, 0.59, 0.874, 0.883, 0.724, 0.34, 0.213, 0.406, 0.978, 0.969, 0.692, 0.205, 0.232, 0.034, 0.27, 0.256, 0.429, 0.721, 0.177, 0.619, 0.71, 0.834, 0.108, 0.689, 0.981, 0.659, 0.254, 0.974, 0.69, 0.23, 0.793, 0.765, 0.012, 0.333, 0.33, 0.713, 0.434, 0.556, 0.683, 0.393, 0.115, 0.616, 0.33, 0.515, 0.286, 0.588, 0.622, 0.505, 0.581, 0.683, 0.737, 0.049, 0.479, 0.344, 0.496, 0.048, 0.773, 0.844, 0.428, 0.901, 0.526, 0.412, 0.013, 0.311, 0.527, 0.74, 0.855, 0.293, 0.519, 0.551, 0.635, 0.732, 0.065, 0.065]
global q = [0.804, 0.618, 0.653, 0.878, 0.248, 0.712, 0.366, 0.903, 0.789, 0.876, 0.942, 0.835, 0.886, 0.973, 0.897, 0.951, 0.975, 0.405, 0.897, 0.522, 0.768, 0.807, 0.824, 0.922, 0.909, 0.527, 0.822, 0.961, 0.626, 0.839, 0.886, 0.89, 0.893, 0.359, 0.899, 0.865, 0.508, 0.947, 0.989, 0.945, 0.935, 0.894, 0.976, 0.265, 0.988, 0.929, 0.883, 0.702, 0.772, 0.261, 0.663, 0.798, 0.978, 0.765, 0.526, 0.893, 0.879, 0.395, 0.988, 0.308, 0.689, 0.757, 0.999, 0.845, 0.9, 0.409, 0.992, 0.812, 0.997, 0.992, 0.698, 0.991, 0.793, 0.795, 0.971, 0.864, 0.799, 0.946, 0.692, 0.599, 0.69, 0.287, 0.874, 0.451, 0.848, 0.896, 0.726, 0.829, 0.923, 0.606, 0.744, 0.894, 0.406, 0.999, 0.996, 0.461, 0.581, 0.506, 0.928, 0.754, 0.964, 0.634, 0.971, 0.649, 0.937, 0.677, 0.955, 0.902, 0.86, 0.979, 0.977, 0.713, 0.523, 0.926, 0.48, 0.898, 0.966, 0.963, 0.919, 0.074, 0.962, 0.861, 0.873, 0.749, 0.85, 0.655, 0.823, 0.847, 0.652, 0.714, 0.89, 0.859, 0.876, 0.958, 0.815, 0.707, 0.452, 0.62, 0.821, 0.948, 0.613, 0.578, 0.939, 0.894, 0.599, 0.985, 0.994, 0.673, 0.54, 0.764, 0.98, 0.791, 0.996, 0.876, 0.299, 0.294, 0.815, 0.96, 0.461, 0.288, 0.959, 0.976, 0.906, 0.903, 0.385, 0.324, 0.849, 0.997, 0.99, 0.697, 0.23, 0.592, 0.814, 0.584, 0.63, 0.914, 0.96, 0.625, 0.726, 0.858, 0.847, 0.348, 0.863, 0.995, 0.988, 0.346, 0.975, 0.889, 0.868, 0.971, 0.997, 0.225, 0.5, 0.965, 0.831, 0.764, 0.95, 0.969, 0.554, 0.86, 0.756, 0.768, 0.849, 0.472, 0.985, 0.888, 0.632, 0.893, 0.867, 0.75, 0.88, 0.941, 0.557, 0.921, 0.823, 0.836, 0.858, 0.468, 0.903, 0.597, 0.856, 0.994, 0.939, 0.802, 0.925, 0.87, 0.669, 0.929, 0.908, 0.711, 0.865, 0.397, 0.295]
global origin = 1
global destination = 50