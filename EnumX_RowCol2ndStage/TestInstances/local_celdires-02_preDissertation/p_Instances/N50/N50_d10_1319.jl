global arcs = [1 25; 1 39; 2 8; 2 37; 2 43; 3 2; 3 6; 3 16; 3 30; 3 43; 3 48; 3 49; 4 9; 4 10; 4 17; 4 23; 4 26; 4 38; 4 43; 4 44; 5 3; 5 9; 5 13; 5 17; 5 25; 5 35; 6 29; 7 15; 7 20; 7 24; 7 27; 7 29; 7 31; 8 5; 8 12; 8 14; 8 28; 8 33; 9 12; 9 16; 9 34; 9 38; 9 47; 10 3; 10 7; 10 9; 10 14; 10 34; 10 36; 10 43; 11 7; 11 9; 11 21; 11 23; 11 28; 11 30; 11 42; 12 3; 12 19; 12 20; 12 36; 12 41; 12 44; 12 49; 13 32; 14 7; 14 9; 14 19; 14 22; 14 36; 14 39; 14 44; 15 10; 15 23; 15 26; 15 44; 15 48; 16 2; 16 7; 16 9; 16 25; 16 27; 16 36; 16 38; 17 18; 17 29; 17 37; 17 46; 18 4; 18 36; 18 38; 18 39; 18 48; 18 49; 19 10; 19 16; 19 24; 20 6; 21 14; 21 18; 21 19; 21 31; 21 49; 22 3; 22 16; 22 17; 22 19; 22 29; 22 31; 22 34; 22 38; 23 5; 23 10; 23 17; 23 27; 23 31; 23 36; 24 12; 24 21; 24 29; 24 34; 25 18; 25 19; 25 32; 25 38; 25 41; 26 7; 26 50; 27 5; 27 6; 27 8; 27 21; 27 34; 28 34; 28 36; 28 45; 28 46; 29 5; 29 7; 29 14; 30 5; 30 7; 30 17; 30 32; 30 47; 31 8; 31 9; 31 19; 31 40; 31 45; 32 4; 32 22; 32 24; 32 26; 32 28; 32 31; 32 39; 33 11; 33 35; 33 36; 33 46; 33 48; 33 50; 34 6; 34 10; 34 12; 34 21; 34 24; 34 38; 35 5; 35 6; 35 39; 36 7; 36 23; 36 24; 36 27; 36 50; 37 4; 37 22; 38 3; 38 8; 38 13; 38 21; 38 34; 38 37; 38 50; 39 9; 39 17; 39 20; 39 21; 40 5; 40 12; 40 26; 40 33; 41 6; 41 22; 41 30; 42 12; 42 18; 42 35; 42 39; 43 9; 43 10; 43 15; 43 17; 43 25; 43 34; 43 46; 44 13; 44 17; 44 34; 44 42; 45 2; 45 4; 45 21; 45 32; 45 34; 46 5; 46 12; 46 13; 46 15; 46 25; 46 32; 47 2; 47 7; 47 11; 47 41; 48 7; 48 36; 48 43; 48 47; 49 19; 49 23; 49 26]
global d_x = [6.0, 1.0, 6.0, 2.0, 5.0, 10.0, 1.0, 9.0, 10.0, 3.0, 5.0, 4.0, 2.0, 5.0, 7.0, 9.0, 4.0, 3.0, 8.0, 9.0, 7.0, 3.0, 4.0, 3.0, 5.0, 5.0, 4.0, 8.0, 4.0, 8.0, 10.0, 1.0, 1.0, 3.0, 7.0, 2.0, 2.0, 9.0, 10.0, 2.0, 10.0, 2.0, 1.0, 7.0, 1.0, 1.0, 3.0, 4.0, 5.0, 10.0, 1.0, 9.0, 7.0, 4.0, 1.0, 10.0, 5.0, 9.0, 7.0, 10.0, 6.0, 9.0, 9.0, 8.0, 7.0, 3.0, 2.0, 10.0, 4.0, 10.0, 4.0, 5.0, 8.0, 3.0, 6.0, 7.0, 9.0, 2.0, 10.0, 9.0, 4.0, 6.0, 3.0, 7.0, 3.0, 5.0, 1.0, 9.0, 3.0, 2.0, 4.0, 9.0, 10.0, 5.0, 2.0, 2.0, 4.0, 1.0, 3.0, 7.0, 7.0, 9.0, 5.0, 3.0, 3.0, 7.0, 10.0, 10.0, 2.0, 2.0, 4.0, 10.0, 4.0, 5.0, 10.0, 2.0, 7.0, 6.0, 9.0, 1.0, 6.0, 6.0, 5.0, 3.0, 9.0, 5.0, 6.0, 3.0, 7.0, 6.0, 1.0, 8.0, 9.0, 9.0, 7.0, 3.0, 6.0, 10.0, 8.0, 5.0, 3.0, 6.0, 10.0, 5.0, 3.0, 7.0, 9.0, 1.0, 9.0, 6.0, 1.0, 5.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0, 2.0, 9.0, 4.0, 2.0, 10.0, 1.0, 3.0, 1.0, 7.0, 2.0, 10.0, 4.0, 5.0, 4.0, 3.0, 5.0, 2.0, 4.0, 9.0, 4.0, 4.0, 7.0, 3.0, 9.0, 2.0, 7.0, 4.0, 4.0, 2.0, 8.0, 9.0, 9.0, 5.0, 1.0, 9.0, 7.0, 2.0, 7.0, 10.0, 9.0, 1.0, 2.0, 10.0, 6.0, 2.0, 5.0, 7.0, 1.0, 8.0, 5.0, 10.0, 2.0, 8.0, 4.0, 2.0, 1.0, 5.0, 6.0, 2.0, 1.0, 7.0, 10.0, 1.0, 1.0, 4.0, 3.0, 3.0, 6.0, 7.0, 3.0, 10.0, 3.0, 2.0, 7.0, 2.0, 4.0]
global b_x = 5
global d_y = [5.0, 6.0, 5.0, 2.0, 2.0, 1.0, 7.0, 2.0, 10.0, 6.0, 3.0, 2.0, 4.0, 9.0, 8.0, 9.0, 4.0, 2.0, 7.0, 10.0, 10.0, 6.0, 6.0, 7.0, 3.0, 3.0, 1.0, 1.0, 8.0, 1.0, 1.0, 6.0, 1.0, 4.0, 4.0, 5.0, 1.0, 7.0, 5.0, 9.0, 5.0, 10.0, 9.0, 5.0, 10.0, 1.0, 1.0, 6.0, 10.0, 6.0, 4.0, 10.0, 4.0, 6.0, 4.0, 9.0, 6.0, 9.0, 8.0, 6.0, 8.0, 2.0, 6.0, 1.0, 6.0, 9.0, 8.0, 7.0, 1.0, 6.0, 3.0, 9.0, 2.0, 6.0, 2.0, 4.0, 8.0, 6.0, 1.0, 8.0, 4.0, 7.0, 7.0, 6.0, 4.0, 5.0, 10.0, 6.0, 5.0, 7.0, 3.0, 8.0, 6.0, 5.0, 4.0, 4.0, 3.0, 2.0, 9.0, 3.0, 3.0, 9.0, 6.0, 2.0, 6.0, 9.0, 9.0, 10.0, 7.0, 2.0, 9.0, 8.0, 3.0, 6.0, 3.0, 6.0, 8.0, 1.0, 7.0, 2.0, 7.0, 7.0, 5.0, 5.0, 8.0, 9.0, 7.0, 10.0, 9.0, 4.0, 10.0, 3.0, 5.0, 5.0, 6.0, 8.0, 8.0, 4.0, 8.0, 5.0, 5.0, 9.0, 3.0, 4.0, 2.0, 8.0, 3.0, 7.0, 8.0, 9.0, 1.0, 9.0, 3.0, 4.0, 2.0, 4.0, 5.0, 3.0, 1.0, 9.0, 8.0, 4.0, 3.0, 3.0, 5.0, 3.0, 2.0, 8.0, 9.0, 3.0, 4.0, 3.0, 4.0, 9.0, 5.0, 4.0, 7.0, 4.0, 9.0, 9.0, 9.0, 9.0, 8.0, 2.0, 7.0, 1.0, 6.0, 4.0, 1.0, 6.0, 2.0, 9.0, 4.0, 1.0, 7.0, 7.0, 2.0, 10.0, 8.0, 3.0, 1.0, 2.0, 5.0, 7.0, 3.0, 10.0, 8.0, 6.0, 7.0, 6.0, 3.0, 1.0, 6.0, 1.0, 1.0, 2.0, 5.0, 6.0, 1.0, 4.0, 8.0, 4.0, 3.0, 8.0, 9.0, 9.0, 7.0, 1.0, 5.0, 8.0, 10.0, 4.0, 10.0, 2.0]
global b_y = 10
global p = [0.449, 0.492, 0.879, 0.328, 0.912, 0.42, 0.08, 0.191, 0.126, 0.429, 0.527, 0.896, 0.413, 0.765, 0.176, 0.67, 0.586, 0.832, 0.985, 0.362, 0.52, 0.594, 0.375, 0.799, 0.597, 0.057, 0.432, 0.885, 0.388, 0.772, 0.389, 0.394, 0.855, 0.077, 0.761, 0.063, 0.6, 0.218, 0.634, 0.007, 0.262, 0.217, 0.763, 0.75, 0.249, 0.469, 0.7, 0.936, 0.437, 0.867, 0.44, 0.55, 0.558, 0.187, 0.08, 0.255, 0.605, 0.551, 0.709, 0.5, 0.603, 0.775, 0.78, 0.753, 0.156, 0.031, 0.578, 0.641, 0.951, 0.201, 0.343, 0.681, 0.539, 0.17, 0.551, 0.265, 0.569, 0.716, 0.032, 0.256, 0.106, 0.519, 0.462, 0.009, 0.461, 0.527, 0.471, 0.691, 0.268, 0.57, 0.072, 0.823, 0.355, 0.535, 0.212, 0.405, 0.293, 0.863, 0.534, 0.733, 0.407, 0.761, 0.461, 0.194, 0.621, 0.303, 0.413, 0.664, 0.555, 0.481, 0.653, 0.42, 0.06, 0.279, 0.786, 0.435, 0.386, 0.652, 0.813, 0.444, 0.968, 0.722, 0.519, 0.469, 0.668, 0.976, 0.207, 0.022, 0.965, 0.129, 0.567, 0.552, 0.386, 0.453, 0.635, 0.409, 0.944, 0.868, 0.846, 0.703, 0.349, 0.012, 0.375, 0.834, 0.442, 0.426, 0.056, 0.51, 0.583, 0.8, 0.697, 0.912, 0.852, 0.833, 0.224, 0.806, 0.519, 0.592, 0.297, 0.357, 0.895, 0.195, 0.048, 0.67, 0.828, 0.65, 0.339, 0.32, 0.898, 0.407, 0.612, 0.768, 0.653, 0.645, 0.077, 0.466, 0.276, 0.306, 0.121, 0.306, 0.367, 0.545, 0.52, 0.606, 0.075, 0.157, 0.309, 0.45, 0.447, 0.227, 0.564, 0.747, 0.697, 0.213, 0.787, 0.108, 0.335, 0.139, 0.1, 0.575, 0.857, 0.636, 0.733, 0.787, 0.617, 0.649, 0.598, 0.51, 0.146, 0.927, 0.453, 0.352, 0.781, 0.609, 0.527, 0.736, 0.107, 0.6, 0.349, 0.502, 0.471, 0.355, 0.381, 0.009, 0.127, 0.002, 0.003, 0.073, 0.029, 0.15, 0.318, 0.773, 0.92, 0.768]
global q = [0.454, 0.778, 0.975, 0.879, 0.923, 0.963, 0.198, 0.306, 0.417, 0.809, 0.711, 0.972, 0.806, 0.948, 0.709, 0.789, 0.806, 0.888, 0.992, 0.811, 0.704, 0.662, 0.727, 0.998, 0.627, 0.356, 0.563, 0.897, 0.779, 0.963, 0.446, 0.865, 0.883, 0.16, 0.879, 0.092, 0.859, 0.548, 0.992, 0.4, 0.664, 0.358, 0.968, 0.943, 0.937, 0.981, 0.721, 0.936, 0.848, 0.928, 0.543, 0.862, 0.653, 0.783, 0.828, 0.422, 0.85, 0.879, 0.892, 0.944, 0.794, 0.778, 0.934, 0.873, 0.396, 0.106, 0.876, 0.99, 0.998, 0.921, 0.794, 0.935, 0.672, 0.521, 0.755, 0.899, 0.857, 0.857, 0.549, 0.631, 0.325, 0.544, 0.769, 0.993, 0.947, 0.565, 0.781, 0.884, 0.373, 0.892, 0.78, 0.854, 0.553, 0.78, 0.526, 0.488, 0.674, 0.938, 0.743, 0.974, 0.973, 0.794, 0.636, 0.433, 0.897, 0.506, 0.72, 0.956, 0.649, 0.816, 0.733, 0.972, 0.927, 0.549, 0.814, 0.968, 0.517, 0.659, 0.813, 0.929, 0.984, 0.997, 0.945, 0.793, 0.89, 0.982, 0.31, 0.159, 0.984, 0.924, 0.78, 0.635, 0.942, 0.644, 0.942, 0.644, 0.997, 0.984, 0.863, 0.886, 0.743, 0.23, 0.688, 0.981, 0.777, 0.65, 0.334, 0.587, 0.615, 0.885, 0.777, 0.995, 0.87, 0.963, 0.41, 0.846, 0.998, 0.753, 0.474, 0.816, 0.927, 0.421, 0.259, 0.9, 0.967, 0.788, 0.816, 0.704, 0.925, 0.414, 0.744, 0.871, 0.999, 0.991, 0.436, 0.602, 0.56, 0.856, 0.999, 0.354, 0.895, 0.724, 0.561, 0.836, 0.155, 0.571, 0.331, 0.77, 0.969, 0.42, 0.945, 0.783, 0.963, 0.772, 0.919, 0.153, 0.502, 0.847, 0.236, 0.822, 0.99, 0.82, 0.959, 0.79, 0.642, 0.94, 0.623, 0.78, 0.491, 0.928, 0.581, 0.911, 0.967, 0.645, 0.764, 0.897, 0.892, 0.97, 0.941, 0.83, 0.74, 0.368, 0.587, 0.055, 0.302, 0.889, 0.914, 0.932, 0.074, 0.886, 0.359, 0.837, 0.937, 0.816]
global origin = 1
global destination = 50