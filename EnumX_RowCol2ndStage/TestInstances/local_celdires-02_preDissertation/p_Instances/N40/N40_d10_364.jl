global arcs = [1 4; 1 9; 1 29; 2 5; 2 9; 2 10; 2 11; 2 22; 2 26; 3 9; 3 20; 3 21; 3 23; 4 3; 4 12; 4 17; 4 24; 4 34; 5 2; 5 11; 5 31; 5 32; 6 7; 6 20; 7 9; 7 16; 7 40; 8 13; 8 29; 8 32; 9 8; 9 39; 10 7; 10 11; 10 17; 10 21; 10 33; 10 35; 10 38; 11 22; 11 23; 11 36; 12 2; 12 11; 12 25; 12 30; 12 35; 13 4; 13 5; 13 26; 13 34; 13 40; 14 8; 14 13; 14 17; 14 24; 14 25; 14 27; 14 39; 15 4; 15 8; 15 11; 15 28; 15 33; 15 39; 16 7; 16 13; 16 21; 17 32; 18 22; 19 8; 19 10; 19 18; 19 40; 20 18; 20 27; 20 34; 21 3; 21 13; 21 27; 21 30; 21 33; 22 4; 22 27; 22 33; 22 37; 23 5; 23 13; 23 29; 23 37; 24 11; 24 15; 24 28; 25 18; 25 29; 25 38; 25 40; 26 3; 26 6; 26 25; 26 33; 26 35; 27 23; 27 25; 27 26; 27 32; 28 4; 28 12; 28 19; 28 25; 28 29; 28 37; 28 38; 29 11; 29 19; 29 21; 29 22; 29 23; 29 39; 30 10; 30 19; 30 20; 31 37; 32 2; 32 15; 32 29; 32 35; 33 2; 33 8; 33 12; 33 13; 33 14; 33 17; 34 5; 34 13; 34 15; 34 20; 34 29; 34 36; 35 22; 35 26; 35 39; 36 4; 36 34; 36 40; 37 8; 37 30; 38 2; 38 27; 38 35; 39 4; 39 12]
global d_x = [2.0, 4.0, 6.0, 2.0, 4.0, 1.0, 5.0, 5.0, 7.0, 1.0, 7.0, 9.0, 5.0, 7.0, 10.0, 6.0, 4.0, 7.0, 5.0, 5.0, 8.0, 6.0, 9.0, 2.0, 10.0, 2.0, 7.0, 9.0, 3.0, 6.0, 9.0, 9.0, 9.0, 2.0, 1.0, 10.0, 4.0, 5.0, 5.0, 1.0, 2.0, 6.0, 5.0, 3.0, 7.0, 2.0, 5.0, 1.0, 2.0, 9.0, 6.0, 8.0, 9.0, 10.0, 4.0, 8.0, 3.0, 7.0, 1.0, 9.0, 2.0, 2.0, 7.0, 10.0, 5.0, 7.0, 7.0, 2.0, 8.0, 8.0, 9.0, 10.0, 4.0, 3.0, 4.0, 1.0, 1.0, 4.0, 7.0, 1.0, 10.0, 3.0, 2.0, 4.0, 10.0, 2.0, 5.0, 8.0, 1.0, 3.0, 10.0, 2.0, 7.0, 8.0, 10.0, 10.0, 6.0, 10.0, 5.0, 1.0, 2.0, 6.0, 2.0, 10.0, 1.0, 4.0, 8.0, 5.0, 8.0, 9.0, 10.0, 8.0, 6.0, 5.0, 7.0, 3.0, 5.0, 9.0, 5.0, 8.0, 3.0, 4.0, 7.0, 7.0, 4.0, 1.0, 8.0, 9.0, 6.0, 2.0, 3.0, 9.0, 2.0, 2.0, 8.0, 7.0, 8.0, 10.0, 2.0, 4.0, 10.0, 1.0, 10.0, 7.0, 4.0, 7.0, 2.0, 1.0, 3.0, 1.0, 2.0, 9.0]
global b_x = 5
global d_y = [7.0, 5.0, 7.0, 1.0, 9.0, 3.0, 2.0, 3.0, 1.0, 10.0, 8.0, 9.0, 4.0, 7.0, 5.0, 2.0, 7.0, 3.0, 6.0, 4.0, 8.0, 4.0, 5.0, 7.0, 8.0, 9.0, 8.0, 6.0, 10.0, 8.0, 10.0, 3.0, 10.0, 10.0, 10.0, 6.0, 7.0, 6.0, 8.0, 5.0, 6.0, 5.0, 4.0, 1.0, 10.0, 9.0, 6.0, 8.0, 1.0, 9.0, 3.0, 3.0, 1.0, 7.0, 8.0, 4.0, 1.0, 7.0, 5.0, 10.0, 7.0, 9.0, 3.0, 8.0, 4.0, 4.0, 1.0, 7.0, 9.0, 5.0, 2.0, 5.0, 7.0, 7.0, 1.0, 6.0, 9.0, 10.0, 7.0, 10.0, 2.0, 3.0, 8.0, 10.0, 4.0, 6.0, 6.0, 5.0, 2.0, 1.0, 9.0, 4.0, 3.0, 1.0, 7.0, 9.0, 10.0, 3.0, 7.0, 8.0, 10.0, 1.0, 10.0, 5.0, 10.0, 10.0, 3.0, 1.0, 5.0, 4.0, 8.0, 9.0, 9.0, 6.0, 7.0, 7.0, 6.0, 8.0, 10.0, 6.0, 2.0, 6.0, 5.0, 2.0, 10.0, 6.0, 6.0, 9.0, 4.0, 5.0, 2.0, 10.0, 8.0, 1.0, 10.0, 7.0, 6.0, 2.0, 8.0, 2.0, 7.0, 6.0, 3.0, 8.0, 3.0, 6.0, 8.0, 10.0, 9.0, 5.0, 6.0, 1.0]
global b_y = 10
global p = [0.428, 0.824, 0.843, 0.694, 0.258, 0.09, 0.608, 0.048, 0.522, 0.584, 0.758, 0.069, 0.32, 0.463, 0.358, 0.153, 0.425, 0.583, 0.409, 0.022, 0.562, 0.299, 0.019, 0.255, 0.707, 0.48, 0.315, 0.876, 0.746, 0.566, 0.16, 0.231, 0.521, 0.906, 0.444, 0.873, 0.764, 0.65, 0.923, 0.34, 0.438, 0.464, 0.026, 0.812, 0.375, 0.054, 0.686, 0.179, 0.056, 0.075, 0.437, 0.293, 0.875, 0.999, 0.955, 0.52, 0.892, 0.351, 0.744, 0.589, 0.757, 0.781, 0.166, 0.21, 0.095, 0.921, 0.59, 0.841, 0.226, 0.011, 0.188, 0.6, 0.962, 0.177, 0.889, 0.962, 0.856, 0.94, 0.905, 0.42, 0.063, 0.515, 0.611, 0.425, 0.1, 0.725, 0.638, 0.487, 0.341, 0.497, 0.446, 0.663, 0.815, 0.153, 0.584, 0.112, 0.366, 0.297, 0.537, 0.051, 0.245, 0.598, 0.664, 0.785, 0.19, 0.672, 0.177, 0.246, 0.794, 0.526, 0.157, 0.573, 0.225, 0.668, 0.098, 0.4, 0.617, 0.764, 0.022, 0.238, 0.894, 0.383, 0.685, 0.806, 0.895, 0.313, 0.946, 0.211, 0.263, 0.146, 0.892, 0.004, 0.973, 0.617, 0.897, 0.998, 0.376, 0.412, 0.846, 0.293, 0.433, 0.365, 0.54, 0.235, 0.708, 0.415, 0.729, 0.645, 0.416, 0.795, 0.67, 0.312]
global q = [0.573, 0.933, 0.87, 0.939, 0.385, 0.473, 0.79, 0.304, 0.692, 0.722, 0.901, 0.094, 0.393, 0.561, 0.37, 0.167, 0.74, 0.62, 0.958, 0.06, 0.89, 0.544, 0.128, 0.681, 0.869, 0.736, 0.819, 0.926, 0.807, 0.925, 0.875, 0.959, 0.7, 0.948, 0.489, 0.887, 0.785, 0.881, 0.943, 0.452, 0.904, 0.576, 0.454, 0.867, 0.918, 0.904, 0.784, 0.444, 0.116, 0.936, 0.587, 0.908, 0.905, 0.999, 0.974, 0.616, 0.98, 0.857, 0.784, 0.956, 0.78, 0.806, 0.503, 0.58, 0.309, 0.988, 0.607, 0.884, 0.566, 0.418, 0.245, 0.649, 0.989, 0.227, 0.952, 0.965, 0.878, 0.968, 0.972, 0.6, 0.848, 0.706, 0.883, 0.639, 0.451, 0.794, 0.945, 0.872, 0.765, 0.891, 0.881, 0.823, 0.872, 0.39, 0.669, 0.628, 0.979, 0.314, 0.899, 0.86, 0.949, 0.964, 0.965, 0.829, 0.251, 0.75, 0.415, 0.517, 0.957, 0.649, 0.584, 0.981, 0.373, 0.699, 0.323, 0.983, 0.913, 0.92, 0.168, 0.283, 0.935, 0.476, 0.806, 0.976, 0.947, 0.696, 0.997, 0.787, 0.301, 0.634, 0.963, 0.741, 0.981, 0.68, 0.956, 0.999, 0.672, 0.928, 0.862, 0.983, 0.817, 0.682, 0.979, 0.719, 0.972, 0.821, 0.982, 0.87, 0.665, 0.945, 0.975, 0.352]
global origin = 1
global destination = 40