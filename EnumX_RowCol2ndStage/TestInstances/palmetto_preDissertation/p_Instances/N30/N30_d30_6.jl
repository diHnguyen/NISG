global arcs = [1 2; 1 3; 1 9; 1 11; 1 12; 1 14; 1 18; 1 19; 1 20; 1 27; 1 30; 2 5; 2 8; 2 19; 2 20; 2 23; 2 24; 3 2; 3 5; 3 6; 3 7; 3 12; 3 16; 3 18; 3 19; 3 20; 3 21; 3 26; 3 28; 3 29; 4 3; 4 7; 4 11; 4 14; 4 17; 4 20; 4 24; 4 25; 4 29; 5 2; 5 6; 5 14; 5 15; 5 27; 6 20; 6 21; 6 27; 6 28; 7 11; 7 17; 7 23; 7 24; 8 4; 8 13; 8 15; 8 17; 8 19; 8 20; 8 24; 8 29; 8 30; 9 3; 9 4; 9 5; 9 7; 9 10; 9 13; 9 15; 9 18; 9 20; 9 23; 9 24; 10 4; 10 5; 10 8; 10 11; 10 12; 10 13; 10 14; 10 16; 10 19; 10 20; 10 22; 10 26; 10 27; 10 30; 11 2; 11 4; 11 17; 11 21; 11 25; 11 26; 11 29; 11 30; 12 3; 12 4; 12 6; 12 7; 12 11; 12 14; 12 18; 12 22; 12 23; 12 24; 13 12; 13 19; 13 21; 13 27; 13 29; 13 30; 14 4; 14 7; 14 15; 14 18; 14 19; 14 20; 14 29; 15 2; 15 6; 15 7; 15 11; 15 12; 15 17; 15 26; 15 27; 15 28; 15 30; 16 5; 16 8; 16 12; 16 13; 16 14; 16 18; 16 21; 16 23; 16 29; 17 5; 17 6; 17 8; 17 9; 17 11; 17 12; 17 13; 17 15; 17 16; 17 25; 17 26; 17 27; 17 29; 17 30; 18 2; 18 5; 18 9; 18 10; 18 11; 18 12; 18 14; 18 16; 18 17; 18 20; 18 23; 18 25; 18 26; 18 28; 18 29; 19 9; 19 12; 19 13; 19 16; 19 24; 19 28; 19 29; 20 2; 20 9; 20 12; 20 17; 20 23; 20 24; 20 27; 21 11; 21 12; 21 14; 21 16; 21 17; 21 23; 21 24; 21 25; 21 27; 21 30; 22 2; 22 4; 22 5; 22 7; 22 9; 22 15; 22 19; 22 24; 22 29; 23 10; 23 25; 23 26; 23 30; 24 2; 24 3; 24 4; 24 8; 24 9; 24 12; 24 13; 24 17; 24 20; 24 26; 24 29; 25 2; 25 7; 25 8; 25 14; 25 19; 25 23; 25 30; 26 7; 26 13; 26 15; 26 16; 26 18; 26 22; 26 23; 26 24; 26 30; 27 2; 27 3; 27 6; 27 8; 27 10; 27 18; 27 22; 27 30; 28 4; 28 5; 28 7; 28 9; 28 13; 28 16; 28 17; 28 27; 29 18; 29 22; 29 27]
global d_x = [5.0, 4.0, 4.0, 7.0, 8.0, 6.0, 9.0, 2.0, 4.0, 5.0, 4.0, 10.0, 3.0, 7.0, 7.0, 2.0, 5.0, 6.0, 7.0, 8.0, 7.0, 1.0, 3.0, 4.0, 10.0, 4.0, 1.0, 3.0, 8.0, 6.0, 7.0, 3.0, 7.0, 8.0, 2.0, 2.0, 6.0, 8.0, 7.0, 3.0, 9.0, 5.0, 9.0, 6.0, 6.0, 6.0, 9.0, 4.0, 5.0, 3.0, 3.0, 6.0, 2.0, 8.0, 3.0, 1.0, 10.0, 10.0, 2.0, 5.0, 10.0, 9.0, 9.0, 4.0, 6.0, 2.0, 1.0, 9.0, 6.0, 9.0, 2.0, 2.0, 10.0, 7.0, 6.0, 8.0, 8.0, 5.0, 2.0, 5.0, 6.0, 8.0, 2.0, 2.0, 7.0, 10.0, 6.0, 9.0, 10.0, 4.0, 1.0, 9.0, 2.0, 8.0, 2.0, 10.0, 1.0, 1.0, 1.0, 8.0, 10.0, 7.0, 4.0, 9.0, 6.0, 3.0, 6.0, 10.0, 3.0, 8.0, 3.0, 4.0, 3.0, 5.0, 2.0, 9.0, 6.0, 4.0, 4.0, 3.0, 8.0, 4.0, 10.0, 7.0, 3.0, 1.0, 10.0, 8.0, 9.0, 9.0, 8.0, 6.0, 9.0, 1.0, 9.0, 6.0, 7.0, 8.0, 2.0, 3.0, 4.0, 10.0, 3.0, 5.0, 5.0, 10.0, 10.0, 10.0, 7.0, 8.0, 8.0, 2.0, 10.0, 7.0, 4.0, 8.0, 7.0, 8.0, 5.0, 1.0, 5.0, 7.0, 3.0, 1.0, 6.0, 10.0, 10.0, 7.0, 6.0, 3.0, 3.0, 2.0, 4.0, 8.0, 1.0, 5.0, 8.0, 7.0, 1.0, 10.0, 9.0, 3.0, 7.0, 8.0, 4.0, 1.0, 3.0, 1.0, 6.0, 7.0, 6.0, 4.0, 1.0, 1.0, 1.0, 1.0, 8.0, 1.0, 7.0, 6.0, 9.0, 8.0, 10.0, 7.0, 8.0, 8.0, 5.0, 8.0, 5.0, 5.0, 1.0, 6.0, 2.0, 2.0, 7.0, 7.0, 7.0, 8.0, 10.0, 5.0, 6.0, 2.0, 1.0, 9.0, 1.0, 5.0, 1.0, 4.0, 9.0, 6.0, 9.0, 4.0, 4.0, 1.0, 9.0, 6.0, 5.0, 5.0, 3.0, 3.0, 2.0, 8.0, 5.0, 6.0, 5.0, 9.0, 1.0, 1.0]
global b_x = 5
global d_y = [3.0, 2.0, 5.0, 9.0, 10.0, 6.0, 9.0, 5.0, 10.0, 10.0, 9.0, 9.0, 4.0, 3.0, 3.0, 2.0, 8.0, 3.0, 7.0, 4.0, 6.0, 4.0, 9.0, 8.0, 8.0, 7.0, 8.0, 7.0, 9.0, 6.0, 9.0, 9.0, 1.0, 4.0, 2.0, 10.0, 6.0, 7.0, 3.0, 2.0, 3.0, 1.0, 7.0, 1.0, 9.0, 10.0, 6.0, 1.0, 1.0, 2.0, 1.0, 1.0, 5.0, 3.0, 3.0, 7.0, 10.0, 9.0, 9.0, 7.0, 4.0, 8.0, 9.0, 10.0, 9.0, 5.0, 8.0, 5.0, 7.0, 6.0, 4.0, 5.0, 2.0, 10.0, 3.0, 4.0, 10.0, 5.0, 2.0, 10.0, 2.0, 10.0, 8.0, 4.0, 3.0, 10.0, 7.0, 3.0, 4.0, 2.0, 6.0, 4.0, 2.0, 5.0, 10.0, 9.0, 9.0, 7.0, 9.0, 8.0, 4.0, 9.0, 3.0, 3.0, 9.0, 9.0, 7.0, 1.0, 6.0, 5.0, 4.0, 2.0, 2.0, 1.0, 4.0, 6.0, 4.0, 5.0, 1.0, 1.0, 9.0, 9.0, 5.0, 6.0, 8.0, 3.0, 7.0, 8.0, 5.0, 7.0, 6.0, 1.0, 10.0, 2.0, 6.0, 2.0, 2.0, 7.0, 4.0, 2.0, 10.0, 2.0, 1.0, 9.0, 7.0, 2.0, 4.0, 8.0, 1.0, 2.0, 6.0, 5.0, 1.0, 7.0, 5.0, 6.0, 3.0, 8.0, 10.0, 2.0, 5.0, 1.0, 8.0, 2.0, 1.0, 4.0, 7.0, 3.0, 5.0, 5.0, 1.0, 10.0, 10.0, 9.0, 8.0, 1.0, 9.0, 9.0, 6.0, 3.0, 8.0, 4.0, 9.0, 7.0, 4.0, 3.0, 9.0, 10.0, 3.0, 3.0, 7.0, 10.0, 6.0, 2.0, 10.0, 8.0, 7.0, 3.0, 4.0, 1.0, 2.0, 4.0, 1.0, 10.0, 3.0, 7.0, 2.0, 9.0, 8.0, 5.0, 1.0, 5.0, 1.0, 1.0, 5.0, 1.0, 2.0, 5.0, 8.0, 6.0, 9.0, 9.0, 8.0, 9.0, 8.0, 8.0, 4.0, 3.0, 3.0, 7.0, 8.0, 2.0, 10.0, 4.0, 8.0, 1.0, 2.0, 8.0, 6.0, 1.0, 3.0, 2.0, 10.0, 7.0, 4.0, 2.0, 6.0, 3.0]
global b_y = 10
global p = [0.035, 0.669, 0.782, 0.175, 0.671, 0.848, 0.207, 0.913, 0.69, 0.104, 0.415, 0.944, 0.57, 0.691, 0.175, 0.514, 0.977, 0.178, 0.508, 0.053, 0.989, 0.987, 0.482, 0.392, 0.964, 0.201, 0.267, 0.22, 0.374, 0.119, 0.608, 0.71, 0.933, 0.332, 0.889, 0.796, 0.706, 0.4, 0.469, 0.834, 0.108, 0.558, 0.281, 0.307, 0.038, 0.778, 0.997, 0.377, 0.953, 0.495, 0.874, 0.107, 0.229, 0.718, 0.522, 0.433, 0.758, 0.67, 0.8, 0.842, 0.51, 0.593, 0.099, 0.749, 0.246, 0.986, 0.69, 0.576, 0.44, 0.015, 0.567, 0.064, 0.435, 0.764, 0.139, 0.121, 0.913, 0.187, 0.284, 0.679, 0.204, 0.284, 0.631, 0.228, 0.793, 0.767, 0.822, 0.176, 0.103, 0.057, 0.624, 0.164, 0.255, 0.922, 0.942, 0.932, 0.235, 0.433, 0.195, 0.151, 0.844, 0.457, 0.796, 0.336, 0.177, 0.847, 0.348, 0.312, 0.2, 0.859, 0.826, 0.355, 0.929, 0.638, 0.331, 0.602, 0.411, 0.67, 0.466, 0.719, 0.584, 0.505, 0.03, 0.649, 0.153, 0.314, 0.041, 0.169, 0.062, 0.036, 0.594, 0.8, 0.183, 0.075, 0.159, 0.157, 0.051, 0.477, 0.356, 0.526, 0.392, 0.677, 0.899, 0.638, 0.352, 0.787, 0.403, 0.046, 0.926, 0.354, 0.267, 0.05, 0.98, 0.718, 0.61, 0.045, 0.484, 0.706, 0.65, 0.95, 0.962, 0.181, 0.977, 0.881, 0.896, 0.667, 0.737, 0.779, 0.912, 0.01, 0.621, 0.214, 0.301, 0.526, 0.233, 0.895, 0.146, 0.666, 0.668, 0.465, 0.56, 0.886, 0.256, 0.751, 0.218, 0.047, 0.296, 0.838, 0.26, 0.067, 0.481, 0.122, 0.425, 0.507, 0.478, 0.538, 0.275, 0.349, 0.738, 0.664, 0.759, 0.177, 0.392, 0.278, 0.831, 0.215, 0.629, 0.554, 0.499, 0.526, 0.104, 0.414, 0.697, 0.725, 0.818, 0.508, 0.847, 0.052, 0.63, 0.212, 0.811, 0.66, 0.104, 0.765, 0.651, 0.606, 0.216, 0.055, 0.454, 0.824, 0.7, 0.887, 0.615, 0.363, 0.718, 0.797, 0.186, 0.05, 0.204, 0.808, 0.289, 0.17, 0.879, 0.161, 0.784, 0.649, 0.017, 0.635]
global q = [0.567, 0.905, 0.89, 0.295, 0.863, 0.85, 0.979, 0.998, 0.882, 0.656, 0.775, 0.983, 0.647, 0.864, 0.639, 0.654, 0.987, 0.328, 0.692, 0.076, 0.998, 0.999, 0.836, 0.854, 0.964, 0.952, 0.274, 0.541, 0.713, 0.122, 0.939, 0.761, 0.978, 0.743, 0.986, 0.947, 0.822, 0.727, 0.839, 0.855, 0.64, 0.947, 0.516, 0.574, 0.972, 0.825, 0.997, 0.592, 0.98, 0.979, 0.971, 0.437, 0.82, 0.749, 0.702, 0.541, 0.85, 0.678, 0.996, 0.982, 0.904, 0.715, 0.657, 0.898, 0.574, 0.996, 0.864, 0.911, 0.646, 0.801, 0.833, 0.156, 0.742, 0.933, 0.897, 0.416, 0.931, 0.439, 0.733, 0.745, 0.64, 0.879, 0.667, 0.571, 0.98, 0.95, 0.923, 0.939, 0.166, 0.254, 0.93, 0.505, 0.579, 0.973, 0.984, 0.958, 0.517, 0.738, 0.84, 0.56, 0.955, 0.669, 0.901, 0.85, 0.841, 0.935, 0.763, 0.596, 0.712, 0.945, 0.989, 0.983, 0.937, 0.859, 0.59, 0.917, 0.678, 0.775, 0.541, 0.991, 0.656, 0.672, 0.783, 0.728, 0.752, 0.353, 0.261, 0.395, 0.688, 0.817, 0.675, 0.888, 0.922, 0.34, 0.4, 0.427, 0.303, 0.786, 0.689, 0.832, 0.646, 0.804, 0.936, 0.997, 0.94, 0.884, 0.865, 0.213, 0.997, 0.61, 0.444, 0.609, 0.981, 0.878, 0.664, 0.562, 0.921, 0.925, 0.799, 0.999, 0.983, 0.921, 0.986, 0.917, 0.917, 0.865, 0.962, 0.974, 0.961, 0.927, 0.719, 0.942, 0.8, 0.56, 0.504, 0.997, 0.641, 0.809, 0.961, 0.991, 0.648, 0.982, 0.991, 0.999, 0.221, 0.4, 0.96, 0.928, 0.81, 0.765, 0.959, 0.429, 0.516, 0.824, 0.958, 0.832, 0.512, 0.591, 0.758, 0.837, 0.915, 0.781, 0.543, 0.609, 0.909, 0.263, 0.849, 0.644, 0.544, 0.771, 0.645, 0.547, 0.822, 0.964, 0.828, 0.571, 0.93, 0.542, 0.998, 0.618, 0.912, 0.936, 0.618, 0.821, 0.978, 0.85, 0.951, 0.517, 0.991, 0.992, 0.881, 0.969, 0.646, 0.623, 0.856, 0.888, 0.683, 0.322, 0.31, 0.848, 0.736, 0.229, 0.982, 0.541, 0.841, 0.986, 0.415, 0.757]
global origin = 1
global destination = 30