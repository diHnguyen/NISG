global arcs = [1 4; 1 8; 1 12; 1 25; 1 27; 1 35; 1 44; 1 45; 2 4; 2 11; 2 14; 2 22; 2 46; 3 2; 3 4; 3 9; 3 21; 3 24; 3 35; 4 11; 4 16; 4 32; 4 44; 5 2; 5 10; 5 34; 5 47; 6 2; 6 10; 6 14; 6 30; 7 5; 7 6; 7 25; 7 32; 7 45; 8 16; 8 22; 8 45; 8 48; 9 6; 9 11; 9 21; 9 28; 9 42; 10 7; 10 21; 10 30; 11 15; 11 28; 11 43; 12 32; 12 41; 13 7; 13 12; 13 23; 13 24; 14 6; 14 12; 14 18; 14 30; 15 6; 15 22; 15 26; 15 29; 16 11; 16 20; 16 22; 16 39; 17 27; 17 46; 17 48; 18 5; 18 17; 18 31; 18 35; 18 39; 18 41; 19 14; 19 47; 20 5; 20 22; 20 46; 21 2; 21 23; 21 36; 21 39; 22 13; 22 23; 22 24; 22 26; 22 39; 22 47; 23 7; 23 10; 23 21; 23 30; 23 31; 23 33; 24 11; 24 35; 24 36; 24 50; 25 4; 25 7; 25 8; 25 10; 25 36; 25 41; 26 25; 26 33; 26 46; 27 2; 27 17; 27 30; 27 40; 28 2; 28 3; 28 11; 28 16; 28 26; 28 34; 29 4; 29 6; 29 10; 29 16; 29 27; 29 30; 29 38; 29 43; 30 23; 30 36; 30 38; 31 11; 31 12; 31 14; 31 16; 31 25; 31 43; 32 36; 32 39; 33 4; 33 15; 33 21; 33 38; 33 42; 34 24; 34 26; 34 39; 34 42; 35 3; 35 6; 35 16; 35 38; 35 42; 36 14; 36 16; 36 17; 36 18; 36 19; 36 22; 36 39; 36 45; 37 15; 37 29; 37 35; 37 42; 37 48; 37 49; 38 14; 38 15; 39 22; 39 23; 39 34; 39 40; 39 41; 40 27; 40 47; 41 33; 41 48; 42 15; 42 16; 42 23; 42 26; 42 35; 43 4; 44 2; 44 22; 44 31; 44 35; 44 42; 44 50; 45 12; 45 14; 45 15; 45 19; 45 27; 45 31; 46 3; 46 27; 46 42; 46 44; 46 47; 47 19; 47 21; 47 37; 47 49; 48 6; 48 10; 48 27; 49 12; 49 17; 49 35]
global d_x = [3.0, 6.0, 5.0, 2.0, 1.0, 7.0, 9.0, 9.0, 1.0, 9.0, 10.0, 3.0, 3.0, 1.0, 4.0, 4.0, 10.0, 8.0, 1.0, 3.0, 10.0, 10.0, 4.0, 5.0, 10.0, 1.0, 1.0, 7.0, 3.0, 7.0, 9.0, 2.0, 10.0, 1.0, 7.0, 1.0, 1.0, 7.0, 5.0, 3.0, 5.0, 6.0, 8.0, 7.0, 8.0, 5.0, 1.0, 6.0, 2.0, 3.0, 1.0, 1.0, 1.0, 7.0, 9.0, 5.0, 9.0, 7.0, 2.0, 5.0, 6.0, 10.0, 10.0, 10.0, 9.0, 9.0, 8.0, 7.0, 5.0, 7.0, 4.0, 2.0, 6.0, 9.0, 9.0, 3.0, 2.0, 5.0, 8.0, 2.0, 8.0, 2.0, 10.0, 8.0, 1.0, 9.0, 1.0, 9.0, 7.0, 4.0, 6.0, 2.0, 2.0, 9.0, 10.0, 9.0, 4.0, 5.0, 5.0, 6.0, 1.0, 1.0, 3.0, 10.0, 3.0, 5.0, 7.0, 9.0, 3.0, 9.0, 5.0, 6.0, 3.0, 2.0, 9.0, 7.0, 10.0, 6.0, 1.0, 9.0, 1.0, 5.0, 6.0, 8.0, 7.0, 2.0, 7.0, 4.0, 9.0, 10.0, 3.0, 3.0, 4.0, 1.0, 6.0, 9.0, 8.0, 4.0, 4.0, 10.0, 10.0, 10.0, 6.0, 5.0, 7.0, 6.0, 5.0, 10.0, 9.0, 9.0, 10.0, 9.0, 8.0, 5.0, 3.0, 10.0, 6.0, 10.0, 1.0, 7.0, 3.0, 10.0, 5.0, 8.0, 9.0, 1.0, 1.0, 8.0, 2.0, 10.0, 1.0, 9.0, 1.0, 6.0, 3.0, 5.0, 10.0, 1.0, 10.0, 3.0, 10.0, 8.0, 8.0, 9.0, 10.0, 5.0, 3.0, 9.0, 5.0, 4.0, 6.0, 5.0, 8.0, 10.0, 10.0, 1.0, 8.0, 8.0, 3.0, 10.0, 4.0, 2.0, 10.0, 9.0, 1.0, 9.0, 9.0, 8.0, 6.0, 9.0, 6.0, 1.0, 3.0]
global b_x = 5
global d_y = [9.0, 3.0, 5.0, 6.0, 3.0, 5.0, 8.0, 4.0, 6.0, 6.0, 3.0, 2.0, 3.0, 3.0, 9.0, 2.0, 4.0, 2.0, 9.0, 8.0, 7.0, 8.0, 10.0, 2.0, 7.0, 7.0, 10.0, 6.0, 7.0, 7.0, 6.0, 5.0, 8.0, 2.0, 4.0, 2.0, 4.0, 4.0, 1.0, 7.0, 1.0, 9.0, 3.0, 2.0, 6.0, 10.0, 5.0, 4.0, 2.0, 2.0, 10.0, 9.0, 2.0, 2.0, 4.0, 6.0, 4.0, 2.0, 6.0, 2.0, 8.0, 8.0, 10.0, 9.0, 8.0, 2.0, 10.0, 4.0, 2.0, 6.0, 8.0, 3.0, 2.0, 10.0, 4.0, 7.0, 5.0, 3.0, 7.0, 7.0, 10.0, 3.0, 6.0, 10.0, 3.0, 5.0, 6.0, 6.0, 4.0, 10.0, 3.0, 5.0, 10.0, 6.0, 1.0, 2.0, 10.0, 4.0, 4.0, 10.0, 10.0, 7.0, 5.0, 8.0, 1.0, 4.0, 1.0, 4.0, 2.0, 5.0, 6.0, 8.0, 6.0, 7.0, 2.0, 1.0, 10.0, 4.0, 6.0, 7.0, 10.0, 6.0, 4.0, 3.0, 4.0, 2.0, 2.0, 9.0, 8.0, 2.0, 2.0, 9.0, 2.0, 4.0, 3.0, 2.0, 3.0, 6.0, 1.0, 9.0, 4.0, 9.0, 4.0, 3.0, 3.0, 2.0, 6.0, 8.0, 9.0, 1.0, 7.0, 3.0, 10.0, 9.0, 6.0, 3.0, 4.0, 8.0, 1.0, 6.0, 9.0, 5.0, 1.0, 1.0, 8.0, 4.0, 9.0, 10.0, 2.0, 10.0, 3.0, 7.0, 7.0, 2.0, 3.0, 1.0, 7.0, 7.0, 9.0, 7.0, 2.0, 10.0, 5.0, 6.0, 5.0, 9.0, 3.0, 10.0, 3.0, 2.0, 6.0, 3.0, 9.0, 9.0, 4.0, 1.0, 4.0, 3.0, 2.0, 4.0, 10.0, 9.0, 10.0, 1.0, 4.0, 4.0, 1.0, 3.0, 4.0, 2.0, 1.0, 2.0, 9.0]
global b_y = 10
global p = [0.736, 0.398, 0.909, 0.271, 0.453, 0.193, 0.445, 0.808, 0.642, 0.876, 0.847, 0.838, 0.888, 0.1, 0.348, 0.192, 0.651, 0.7, 0.583, 0.791, 0.806, 0.235, 0.285, 0.46, 0.846, 0.688, 0.991, 0.23, 0.711, 0.537, 0.727, 0.318, 0.086, 0.118, 0.779, 0.554, 0.265, 0.893, 0.81, 0.636, 0.756, 0.42, 0.062, 0.106, 0.075, 0.113, 0.949, 0.816, 0.57, 0.863, 0.762, 0.686, 0.908, 0.55, 0.048, 0.844, 0.377, 0.305, 0.354, 0.925, 0.255, 0.022, 0.822, 0.76, 0.833, 0.009, 0.962, 0.993, 0.676, 0.305, 0.425, 0.885, 0.553, 0.635, 0.058, 0.417, 0.539, 0.033, 0.014, 0.453, 0.823, 0.06, 0.389, 0.923, 0.649, 0.682, 0.273, 0.256, 0.438, 0.243, 0.365, 0.887, 0.464, 0.76, 0.869, 0.349, 0.838, 0.01, 0.135, 0.179, 0.659, 0.908, 0.576, 0.054, 0.078, 0.308, 0.231, 0.972, 0.461, 0.538, 0.192, 0.744, 0.39, 0.174, 0.748, 0.757, 0.507, 0.122, 0.606, 0.434, 0.724, 0.75, 0.746, 0.66, 0.969, 0.551, 0.928, 0.011, 0.096, 0.589, 0.421, 0.297, 0.71, 0.538, 0.328, 0.374, 0.963, 0.043, 0.711, 0.7, 0.09, 0.27, 0.222, 0.391, 0.157, 0.673, 0.549, 0.233, 0.098, 0.056, 0.293, 0.362, 0.575, 0.453, 0.308, 0.911, 0.236, 0.984, 0.771, 0.213, 0.634, 0.549, 0.89, 0.718, 0.746, 0.335, 0.95, 0.801, 0.32, 0.16, 0.888, 0.929, 0.61, 0.532, 0.463, 0.131, 0.526, 0.965, 0.492, 0.796, 0.072, 0.925, 0.72, 0.244, 0.075, 0.402, 0.206, 0.991, 0.445, 0.299, 0.206, 0.272, 0.708, 0.234, 0.858, 0.137, 0.002, 0.397, 0.88, 0.614, 0.998, 0.525, 0.719, 0.459, 0.485, 0.038, 0.377, 0.237, 0.635, 0.516, 0.119, 0.824, 0.658]
global q = [0.909, 0.532, 0.977, 0.507, 0.506, 0.546, 0.805, 0.894, 0.874, 0.906, 0.954, 0.864, 0.921, 0.696, 0.508, 0.637, 0.864, 0.804, 0.7, 0.951, 0.879, 0.664, 0.423, 0.71, 0.916, 0.755, 0.991, 0.461, 0.781, 0.904, 0.763, 0.49, 0.257, 0.41, 0.853, 0.75, 0.84, 0.971, 0.823, 0.944, 0.822, 0.571, 0.1, 0.955, 0.164, 0.833, 0.995, 0.961, 0.741, 0.918, 0.913, 0.909, 0.931, 0.612, 0.365, 0.855, 0.729, 0.99, 0.436, 0.966, 0.628, 0.355, 0.969, 0.796, 0.986, 0.372, 0.974, 0.994, 0.903, 0.778, 0.755, 0.947, 0.604, 0.664, 0.833, 0.799, 0.724, 0.232, 0.672, 0.826, 0.91, 0.934, 0.864, 0.949, 0.7, 0.998, 0.658, 0.487, 0.866, 0.92, 0.685, 0.9, 0.784, 0.985, 0.912, 0.95, 0.898, 0.459, 0.631, 0.671, 0.761, 0.945, 0.652, 0.274, 0.624, 0.438, 0.7, 0.988, 0.82, 0.538, 0.307, 0.982, 0.871, 0.581, 0.749, 0.799, 0.955, 0.269, 0.823, 0.55, 0.756, 0.973, 0.988, 0.809, 0.991, 0.968, 0.934, 0.384, 0.768, 0.697, 0.655, 0.822, 0.842, 0.887, 0.471, 0.556, 0.999, 0.818, 0.965, 0.75, 0.929, 0.649, 0.934, 0.991, 0.567, 0.692, 0.944, 0.403, 0.57, 0.913, 0.999, 0.85, 0.85, 0.919, 0.925, 0.98, 0.61, 0.984, 0.881, 0.846, 0.971, 0.921, 0.912, 0.768, 0.907, 0.5, 0.964, 0.864, 0.652, 0.384, 0.959, 0.956, 0.955, 0.586, 0.889, 0.983, 0.862, 0.976, 0.588, 0.904, 0.938, 0.962, 0.975, 0.525, 0.527, 0.46, 0.228, 0.992, 0.576, 0.631, 0.437, 0.879, 0.837, 0.711, 0.881, 0.217, 0.334, 0.81, 0.955, 0.672, 0.999, 0.668, 0.979, 0.527, 0.929, 0.144, 0.806, 0.664, 0.714, 0.832, 0.136, 0.958, 0.9]
global origin = 1
global destination = 50