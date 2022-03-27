global arcs = [1 3; 1 5; 1 7; 1 9; 1 10; 1 14; 2 9; 2 27; 2 36; 2 38; 3 5; 4 21; 5 33; 5 35; 5 39; 6 5; 6 7; 6 10; 6 23; 6 26; 6 27; 6 35; 6 40; 7 3; 7 6; 7 10; 7 16; 7 17; 8 7; 8 11; 8 15; 8 20; 8 32; 8 33; 9 2; 9 3; 9 10; 9 12; 9 13; 9 14; 10 3; 10 4; 10 13; 10 17; 11 13; 11 23; 11 33; 12 4; 12 7; 12 8; 12 11; 12 14; 12 16; 12 34; 13 28; 13 31; 14 4; 14 13; 14 15; 14 26; 14 27; 14 28; 14 34; 15 8; 15 9; 15 22; 15 26; 15 27; 15 28; 15 31; 15 33; 15 34; 16 6; 16 8; 16 10; 16 14; 16 17; 16 33; 16 35; 16 40; 17 3; 17 29; 17 34; 17 39; 18 25; 18 27; 19 18; 20 9; 20 11; 20 31; 21 9; 21 11; 21 20; 21 26; 21 33; 21 37; 22 5; 22 6; 22 15; 22 16; 23 11; 23 28; 23 29; 23 30; 24 14; 24 17; 24 26; 24 36; 25 2; 25 17; 25 27; 26 14; 26 15; 26 17; 26 22; 26 24; 26 28; 26 31; 26 34; 27 16; 28 9; 28 10; 28 21; 29 7; 29 11; 29 33; 29 39; 30 13; 30 18; 30 32; 30 36; 31 15; 31 23; 31 38; 32 3; 32 11; 32 21; 32 25; 33 16; 34 3; 34 4; 34 11; 34 19; 34 31; 34 35; 35 14; 35 16; 35 17; 35 33; 36 19; 37 7; 37 16; 37 34; 37 36; 38 2; 38 3; 38 19; 38 32; 38 33; 39 6; 39 8; 39 11; 39 18; 39 29]
global d_x = [6.0, 7.0, 6.0, 8.0, 8.0, 4.0, 4.0, 2.0, 5.0, 3.0, 3.0, 2.0, 10.0, 4.0, 10.0, 3.0, 9.0, 1.0, 1.0, 7.0, 5.0, 3.0, 6.0, 7.0, 9.0, 2.0, 1.0, 7.0, 9.0, 6.0, 6.0, 8.0, 3.0, 8.0, 6.0, 2.0, 7.0, 7.0, 7.0, 8.0, 8.0, 1.0, 1.0, 9.0, 3.0, 7.0, 9.0, 8.0, 3.0, 6.0, 4.0, 9.0, 8.0, 7.0, 9.0, 1.0, 2.0, 5.0, 10.0, 1.0, 8.0, 2.0, 5.0, 5.0, 5.0, 8.0, 3.0, 2.0, 9.0, 4.0, 9.0, 7.0, 8.0, 10.0, 8.0, 10.0, 6.0, 4.0, 7.0, 3.0, 1.0, 4.0, 2.0, 8.0, 1.0, 1.0, 3.0, 5.0, 2.0, 10.0, 9.0, 9.0, 10.0, 4.0, 10.0, 5.0, 5.0, 6.0, 2.0, 5.0, 7.0, 8.0, 10.0, 8.0, 9.0, 1.0, 3.0, 4.0, 1.0, 10.0, 1.0, 4.0, 2.0, 2.0, 9.0, 2.0, 2.0, 5.0, 2.0, 7.0, 5.0, 8.0, 4.0, 2.0, 3.0, 1.0, 2.0, 7.0, 6.0, 5.0, 8.0, 6.0, 10.0, 7.0, 5.0, 2.0, 1.0, 1.0, 2.0, 2.0, 9.0, 6.0, 1.0, 8.0, 1.0, 4.0, 6.0, 4.0, 5.0, 6.0, 2.0, 9.0, 2.0, 4.0, 1.0, 1.0, 5.0, 7.0, 9.0, 5.0, 7.0, 5.0, 5.0, 4.0]
global b_x = 5
global d_y = [10.0, 3.0, 6.0, 3.0, 7.0, 6.0, 8.0, 1.0, 10.0, 6.0, 1.0, 9.0, 10.0, 5.0, 1.0, 4.0, 6.0, 10.0, 6.0, 10.0, 1.0, 7.0, 2.0, 9.0, 2.0, 6.0, 1.0, 7.0, 8.0, 7.0, 10.0, 10.0, 2.0, 8.0, 8.0, 3.0, 9.0, 9.0, 10.0, 5.0, 2.0, 1.0, 9.0, 9.0, 2.0, 8.0, 2.0, 8.0, 10.0, 1.0, 3.0, 10.0, 4.0, 3.0, 5.0, 9.0, 2.0, 10.0, 9.0, 6.0, 2.0, 9.0, 5.0, 3.0, 2.0, 1.0, 8.0, 9.0, 1.0, 1.0, 3.0, 2.0, 9.0, 7.0, 5.0, 7.0, 6.0, 9.0, 2.0, 3.0, 1.0, 10.0, 4.0, 4.0, 6.0, 5.0, 1.0, 10.0, 5.0, 5.0, 2.0, 10.0, 3.0, 2.0, 3.0, 5.0, 3.0, 3.0, 7.0, 2.0, 10.0, 2.0, 6.0, 7.0, 10.0, 2.0, 8.0, 8.0, 9.0, 2.0, 5.0, 2.0, 5.0, 5.0, 9.0, 2.0, 6.0, 4.0, 10.0, 3.0, 2.0, 7.0, 1.0, 4.0, 9.0, 5.0, 1.0, 8.0, 1.0, 8.0, 10.0, 10.0, 7.0, 3.0, 4.0, 4.0, 4.0, 3.0, 6.0, 4.0, 8.0, 4.0, 3.0, 9.0, 2.0, 6.0, 5.0, 2.0, 10.0, 7.0, 9.0, 7.0, 2.0, 4.0, 3.0, 6.0, 5.0, 9.0, 7.0, 4.0, 8.0, 2.0, 6.0, 5.0]
global b_y = 10
global p = [0.151, 0.315, 0.007, 0.922, 0.75, 0.484, 0.587, 0.361, 0.757, 0.06, 0.352, 0.704, 0.963, 0.023, 0.606, 0.119, 0.946, 0.013, 0.984, 0.414, 0.856, 0.684, 0.339, 0.906, 0.787, 0.772, 0.662, 0.456, 0.999, 0.544, 0.497, 0.448, 0.838, 0.884, 0.64, 0.453, 0.93, 0.28, 0.015, 0.544, 0.379, 0.486, 0.419, 0.431, 0.468, 0.515, 0.196, 0.776, 0.685, 0.021, 0.892, 0.067, 0.242, 0.551, 0.069, 0.337, 0.929, 0.813, 0.284, 0.348, 0.864, 0.673, 0.398, 0.905, 0.97, 0.607, 0.542, 0.638, 0.496, 0.849, 0.37, 0.134, 0.614, 0.472, 0.135, 0.671, 0.341, 0.209, 0.278, 0.881, 0.559, 0.435, 0.243, 0.938, 0.813, 0.444, 0.855, 0.165, 0.931, 0.597, 0.855, 0.444, 0.116, 0.4, 0.812, 0.868, 0.977, 0.226, 0.011, 0.812, 0.374, 0.635, 0.758, 0.239, 0.805, 0.675, 0.214, 0.754, 0.789, 0.003, 0.176, 0.641, 0.393, 0.798, 0.198, 0.92, 0.942, 0.265, 0.693, 0.096, 0.789, 0.978, 0.117, 0.636, 0.445, 0.597, 0.497, 0.36, 0.595, 0.865, 0.634, 0.356, 0.179, 0.618, 0.888, 0.32, 0.443, 0.525, 0.468, 0.928, 0.313, 0.384, 0.126, 0.555, 0.579, 0.251, 0.46, 0.995, 0.578, 0.028, 0.567, 0.174, 0.095, 0.796, 0.072, 0.59, 0.42, 0.041, 0.077, 0.327, 0.302, 0.873, 0.922, 0.094]
global q = [0.668, 0.328, 0.794, 0.992, 0.751, 0.876, 0.721, 0.651, 0.985, 0.121, 0.547, 0.802, 0.97, 0.695, 0.916, 0.15, 0.974, 0.278, 0.986, 0.763, 0.91, 0.799, 0.845, 0.97, 0.968, 0.834, 0.702, 0.615, 0.999, 0.815, 0.571, 0.672, 0.855, 0.884, 0.68, 0.485, 0.952, 0.757, 0.508, 0.556, 0.725, 0.885, 0.636, 0.701, 0.497, 0.993, 0.213, 0.782, 0.711, 0.298, 0.932, 0.654, 0.434, 0.573, 0.656, 0.745, 0.932, 0.859, 0.88, 0.475, 0.887, 0.704, 0.779, 0.958, 0.987, 0.873, 0.646, 0.926, 0.62, 0.988, 0.637, 0.319, 0.866, 0.699, 0.4, 0.902, 0.591, 0.398, 0.53, 0.903, 0.979, 0.864, 0.857, 0.976, 0.944, 0.798, 0.986, 0.733, 0.952, 0.736, 0.942, 0.792, 0.708, 0.826, 0.867, 0.917, 0.987, 0.75, 0.362, 0.916, 0.994, 0.85, 0.861, 0.769, 0.968, 0.69, 0.857, 0.807, 0.856, 0.648, 0.864, 0.674, 0.415, 0.849, 0.226, 0.947, 0.951, 0.446, 0.989, 0.24, 0.844, 0.991, 0.875, 0.86, 0.875, 0.84, 0.776, 0.708, 0.765, 0.885, 0.816, 0.488, 0.604, 0.712, 0.987, 0.666, 0.97, 0.985, 0.824, 0.957, 0.795, 0.541, 0.885, 0.637, 0.718, 0.359, 0.887, 0.997, 0.607, 0.298, 0.645, 0.734, 0.532, 0.814, 0.132, 0.806, 0.605, 0.51, 0.708, 0.494, 0.344, 0.947, 0.967, 0.17]
global origin = 1
global destination = 40