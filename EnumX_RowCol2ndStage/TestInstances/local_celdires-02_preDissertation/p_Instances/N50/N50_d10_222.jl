global arcs = [1 9; 1 22; 1 24; 1 38; 1 45; 2 6; 2 16; 2 34; 3 4; 3 15; 3 32; 3 33; 3 39; 4 10; 4 16; 4 20; 4 23; 4 27; 4 32; 4 48; 5 4; 5 14; 5 15; 5 26; 5 27; 5 32; 5 37; 5 42; 5 44; 6 2; 6 4; 6 21; 6 25; 6 28; 6 50; 7 14; 7 45; 7 50; 8 19; 8 23; 8 36; 8 44; 8 48; 9 4; 9 28; 9 36; 10 7; 10 11; 10 12; 10 23; 10 49; 11 28; 12 3; 12 13; 12 33; 12 42; 13 16; 13 17; 13 45; 13 48; 13 49; 14 16; 14 32; 14 40; 14 45; 15 5; 15 12; 15 19; 15 31; 15 32; 15 39; 15 45; 16 23; 16 30; 17 16; 17 33; 17 39; 18 30; 18 40; 18 49; 19 12; 19 20; 19 28; 19 39; 19 43; 20 14; 20 21; 20 34; 20 43; 21 14; 21 27; 21 31; 21 38; 22 7; 22 8; 22 21; 22 24; 22 36; 22 40; 22 48; 23 12; 23 15; 23 25; 23 31; 23 45; 23 47; 24 4; 24 6; 24 14; 24 16; 24 21; 24 29; 24 32; 24 38; 24 41; 24 42; 24 44; 25 3; 25 6; 25 13; 25 14; 25 23; 25 34; 25 37; 25 39; 25 40; 25 47; 25 48; 26 8; 26 10; 26 11; 26 18; 26 25; 26 37; 27 6; 27 26; 27 43; 28 6; 28 10; 28 21; 28 41; 28 46; 29 2; 29 6; 29 25; 30 42; 30 43; 31 17; 31 27; 31 30; 31 39; 32 8; 32 9; 32 11; 32 12; 32 13; 32 14; 32 24; 32 25; 32 45; 33 6; 33 9; 33 12; 33 18; 33 28; 33 29; 34 3; 34 26; 34 32; 34 43; 34 48; 35 20; 35 22; 35 23; 35 27; 35 40; 35 47; 35 50; 36 4; 36 43; 36 44; 36 49; 37 13; 37 20; 37 32; 38 6; 38 24; 38 25; 38 37; 38 49; 39 6; 39 16; 39 20; 39 25; 39 31; 39 42; 40 3; 40 12; 40 24; 40 25; 40 41; 40 44; 41 10; 41 11; 41 16; 41 19; 41 28; 41 30; 41 34; 41 43; 42 2; 42 12; 42 13; 42 39; 42 43; 42 46; 43 5; 43 24; 43 33; 43 39; 44 20; 44 34; 45 24; 45 33; 46 28; 46 31; 46 37; 46 42; 46 44; 47 21; 47 30; 47 39; 47 40; 47 44; 48 4; 49 13; 49 14; 49 17; 49 23; 49 26; 49 29; 49 31; 49 32; 49 34; 49 35; 49 48; 49 50]
global d_x = [6.0, 1.0, 2.0, 8.0, 9.0, 9.0, 9.0, 3.0, 1.0, 2.0, 6.0, 7.0, 8.0, 3.0, 8.0, 4.0, 8.0, 5.0, 8.0, 7.0, 5.0, 3.0, 2.0, 2.0, 10.0, 3.0, 8.0, 5.0, 10.0, 9.0, 7.0, 2.0, 2.0, 10.0, 10.0, 1.0, 2.0, 8.0, 4.0, 8.0, 7.0, 10.0, 1.0, 3.0, 7.0, 8.0, 8.0, 9.0, 3.0, 3.0, 5.0, 10.0, 2.0, 1.0, 5.0, 10.0, 9.0, 9.0, 10.0, 3.0, 4.0, 9.0, 7.0, 10.0, 2.0, 8.0, 10.0, 8.0, 5.0, 5.0, 10.0, 1.0, 4.0, 4.0, 5.0, 10.0, 9.0, 7.0, 6.0, 10.0, 10.0, 5.0, 5.0, 10.0, 1.0, 10.0, 9.0, 2.0, 3.0, 6.0, 1.0, 7.0, 10.0, 8.0, 3.0, 10.0, 8.0, 9.0, 3.0, 4.0, 9.0, 8.0, 10.0, 4.0, 10.0, 7.0, 5.0, 4.0, 8.0, 8.0, 9.0, 10.0, 5.0, 9.0, 2.0, 6.0, 8.0, 7.0, 3.0, 8.0, 7.0, 3.0, 1.0, 1.0, 3.0, 1.0, 6.0, 2.0, 3.0, 9.0, 5.0, 9.0, 6.0, 7.0, 6.0, 9.0, 5.0, 1.0, 7.0, 4.0, 1.0, 2.0, 8.0, 5.0, 2.0, 6.0, 10.0, 9.0, 4.0, 4.0, 10.0, 8.0, 4.0, 2.0, 8.0, 8.0, 6.0, 8.0, 3.0, 7.0, 3.0, 10.0, 5.0, 8.0, 8.0, 2.0, 1.0, 5.0, 2.0, 9.0, 8.0, 6.0, 8.0, 9.0, 9.0, 8.0, 1.0, 6.0, 7.0, 4.0, 1.0, 5.0, 9.0, 8.0, 3.0, 2.0, 3.0, 10.0, 8.0, 9.0, 10.0, 2.0, 4.0, 2.0, 10.0, 2.0, 10.0, 1.0, 3.0, 3.0, 9.0, 1.0, 10.0, 5.0, 4.0, 3.0, 5.0, 1.0, 5.0, 9.0, 7.0, 8.0, 3.0, 4.0, 5.0, 5.0, 7.0, 3.0, 5.0, 6.0, 3.0, 10.0, 4.0, 1.0, 10.0, 5.0, 4.0, 3.0, 9.0, 9.0, 3.0, 6.0, 2.0, 1.0, 2.0, 2.0, 10.0, 10.0, 9.0, 10.0, 7.0, 1.0, 7.0, 5.0, 4.0, 9.0, 9.0]
global b_x = 5
global d_y = [2.0, 2.0, 6.0, 5.0, 7.0, 9.0, 3.0, 10.0, 2.0, 6.0, 7.0, 3.0, 2.0, 8.0, 2.0, 2.0, 10.0, 7.0, 8.0, 5.0, 10.0, 7.0, 9.0, 6.0, 9.0, 3.0, 3.0, 5.0, 6.0, 1.0, 7.0, 10.0, 8.0, 2.0, 7.0, 6.0, 7.0, 4.0, 6.0, 5.0, 7.0, 7.0, 5.0, 7.0, 7.0, 7.0, 1.0, 9.0, 3.0, 4.0, 9.0, 2.0, 5.0, 6.0, 9.0, 5.0, 1.0, 7.0, 2.0, 5.0, 3.0, 2.0, 10.0, 5.0, 5.0, 8.0, 8.0, 7.0, 1.0, 2.0, 4.0, 3.0, 6.0, 9.0, 10.0, 9.0, 8.0, 4.0, 3.0, 8.0, 6.0, 3.0, 2.0, 5.0, 1.0, 6.0, 3.0, 10.0, 9.0, 10.0, 4.0, 5.0, 10.0, 1.0, 9.0, 4.0, 5.0, 10.0, 9.0, 5.0, 6.0, 3.0, 7.0, 8.0, 9.0, 4.0, 8.0, 10.0, 7.0, 7.0, 5.0, 10.0, 6.0, 1.0, 3.0, 3.0, 7.0, 7.0, 6.0, 4.0, 1.0, 2.0, 6.0, 4.0, 8.0, 4.0, 6.0, 9.0, 7.0, 2.0, 7.0, 7.0, 8.0, 10.0, 6.0, 1.0, 9.0, 10.0, 2.0, 6.0, 8.0, 5.0, 5.0, 1.0, 5.0, 7.0, 9.0, 9.0, 6.0, 6.0, 9.0, 9.0, 6.0, 6.0, 6.0, 1.0, 6.0, 3.0, 8.0, 3.0, 3.0, 8.0, 5.0, 8.0, 6.0, 2.0, 4.0, 1.0, 8.0, 5.0, 3.0, 2.0, 7.0, 1.0, 2.0, 6.0, 8.0, 1.0, 5.0, 6.0, 5.0, 5.0, 6.0, 7.0, 7.0, 7.0, 5.0, 1.0, 10.0, 9.0, 4.0, 6.0, 2.0, 2.0, 5.0, 10.0, 10.0, 4.0, 2.0, 10.0, 10.0, 3.0, 7.0, 4.0, 7.0, 9.0, 1.0, 4.0, 7.0, 6.0, 5.0, 1.0, 2.0, 8.0, 5.0, 9.0, 8.0, 5.0, 5.0, 6.0, 8.0, 3.0, 1.0, 2.0, 6.0, 4.0, 3.0, 2.0, 5.0, 10.0, 7.0, 5.0, 6.0, 1.0, 1.0, 7.0, 6.0, 3.0, 5.0, 10.0, 5.0, 8.0, 5.0, 6.0, 5.0, 4.0, 4.0]
global b_y = 10
global p = [0.908, 0.224, 0.038, 0.54, 0.409, 0.761, 0.743, 0.594, 0.74, 0.958, 0.113, 0.743, 0.651, 0.937, 0.715, 0.889, 0.097, 0.744, 0.454, 0.333, 0.866, 0.957, 0.084, 0.041, 0.087, 0.425, 0.567, 0.799, 0.659, 0.693, 0.649, 0.627, 0.117, 0.025, 0.981, 0.285, 0.434, 0.152, 0.589, 0.468, 0.667, 0.738, 0.128, 0.776, 0.835, 0.117, 0.232, 0.834, 0.279, 0.645, 0.006, 0.94, 0.245, 0.552, 0.708, 0.676, 0.084, 0.797, 0.909, 0.597, 0.557, 0.646, 0.98, 0.621, 0.323, 0.122, 0.197, 0.932, 0.404, 0.183, 0.44, 0.887, 0.34, 0.499, 0.703, 0.635, 0.172, 0.373, 0.313, 0.453, 0.065, 0.854, 0.235, 0.449, 0.747, 0.283, 0.553, 0.586, 0.666, 0.385, 0.674, 0.379, 0.747, 0.404, 0.928, 0.608, 0.517, 0.272, 0.439, 0.119, 0.649, 0.888, 0.072, 0.407, 0.24, 0.021, 0.236, 0.162, 0.143, 0.679, 0.779, 0.823, 0.895, 0.679, 0.71, 0.468, 0.752, 0.554, 0.563, 0.077, 0.226, 0.344, 0.66, 0.208, 0.776, 0.418, 0.38, 0.768, 0.732, 0.627, 0.542, 0.787, 0.55, 0.599, 0.236, 0.116, 0.703, 0.254, 0.356, 0.209, 0.84, 0.355, 0.199, 0.618, 0.285, 0.873, 0.232, 0.599, 0.27, 0.452, 0.61, 0.763, 0.51, 0.531, 0.117, 0.854, 0.711, 0.211, 0.476, 0.236, 0.192, 0.684, 0.827, 0.221, 0.944, 0.082, 0.484, 0.072, 0.144, 0.425, 0.581, 0.351, 0.388, 0.833, 0.408, 0.443, 0.939, 0.009, 0.662, 0.972, 0.587, 0.743, 0.015, 0.193, 0.928, 0.989, 0.02, 0.06, 0.844, 0.11, 0.473, 0.963, 0.404, 0.898, 0.794, 0.068, 0.619, 0.842, 0.913, 0.652, 0.438, 0.825, 0.457, 0.05, 0.008, 0.579, 0.549, 0.198, 0.332, 0.712, 0.191, 0.431, 0.42, 0.336, 0.944, 0.607, 0.359, 0.827, 0.482, 0.497, 0.927, 0.071, 0.663, 0.117, 0.363, 0.8, 0.948, 0.324, 0.86, 0.065, 0.575, 0.696, 0.482, 0.17, 0.1, 0.126, 0.229, 0.024, 0.525, 0.561, 0.36, 0.968, 0.107, 0.665, 0.911, 0.921, 0.817]
global q = [0.909, 0.671, 0.37, 0.82, 0.753, 0.914, 0.828, 0.871, 0.997, 0.981, 0.614, 0.964, 0.656, 0.963, 0.927, 0.98, 0.159, 0.783, 0.808, 0.748, 0.878, 0.975, 0.797, 0.745, 0.826, 0.858, 0.639, 0.924, 0.817, 0.837, 0.97, 0.721, 0.223, 0.902, 0.993, 0.823, 0.538, 0.982, 0.631, 0.685, 0.798, 0.937, 0.739, 0.881, 0.963, 0.855, 0.384, 0.881, 0.492, 0.995, 0.021, 0.994, 0.613, 0.852, 0.967, 0.992, 0.946, 0.874, 0.987, 0.953, 0.574, 0.879, 0.984, 0.929, 0.5, 0.252, 0.662, 0.937, 0.414, 0.378, 0.887, 0.921, 0.909, 0.539, 0.844, 0.834, 0.511, 0.511, 0.726, 0.908, 0.782, 0.955, 0.967, 0.794, 0.948, 0.585, 0.676, 0.978, 0.73, 0.888, 0.903, 0.656, 0.905, 0.431, 0.971, 0.669, 0.709, 0.631, 0.74, 0.173, 0.943, 0.991, 0.662, 0.966, 0.538, 0.701, 0.693, 0.483, 0.666, 0.982, 0.929, 0.92, 0.945, 0.711, 0.882, 0.69, 0.902, 0.707, 0.704, 0.502, 0.234, 0.585, 0.98, 0.848, 0.913, 0.712, 0.693, 0.814, 0.825, 0.794, 0.573, 0.986, 0.847, 0.883, 0.56, 0.117, 0.952, 0.517, 0.815, 0.807, 0.842, 0.642, 0.514, 0.928, 0.769, 0.99, 0.682, 0.656, 0.746, 0.491, 0.848, 0.965, 0.624, 0.669, 0.196, 0.888, 0.809, 0.694, 0.575, 0.344, 0.254, 0.804, 0.93, 0.536, 0.989, 0.406, 0.686, 0.917, 0.9, 0.583, 0.883, 0.465, 0.513, 0.834, 0.496, 0.523, 0.966, 0.158, 0.986, 0.981, 0.612, 0.775, 0.06, 0.609, 0.937, 0.997, 0.371, 0.621, 0.869, 0.956, 0.748, 0.993, 0.57, 0.899, 0.897, 0.888, 0.794, 0.981, 0.931, 0.922, 0.495, 0.887, 0.924, 0.734, 0.401, 0.949, 0.856, 0.33, 0.773, 0.755, 0.364, 0.496, 0.481, 0.722, 0.962, 0.746, 0.947, 0.827, 0.589, 0.522, 0.946, 0.122, 0.757, 0.305, 0.709, 0.815, 0.976, 0.401, 0.946, 0.104, 0.707, 0.963, 0.566, 0.98, 0.247, 0.59, 0.43, 0.537, 0.862, 0.958, 0.444, 0.978, 0.672, 0.725, 0.916, 0.937, 0.997]
global origin = 1
global destination = 50