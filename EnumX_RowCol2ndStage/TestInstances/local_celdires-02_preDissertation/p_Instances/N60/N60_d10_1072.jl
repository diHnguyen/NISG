global arcs = [1 12; 1 22; 1 42; 1 49; 2 38; 2 59; 3 2; 3 8; 3 25; 3 29; 3 40; 3 60; 4 2; 4 15; 4 44; 4 55; 5 10; 5 13; 5 18; 5 28; 5 48; 5 52; 5 53; 5 55; 5 57; 6 8; 6 17; 6 19; 6 22; 6 25; 6 36; 7 5; 7 12; 7 13; 7 22; 7 30; 7 34; 7 57; 7 58; 8 20; 8 30; 8 40; 8 48; 9 2; 9 18; 9 23; 9 35; 10 20; 10 26; 10 32; 10 36; 10 42; 11 2; 11 5; 11 14; 11 25; 11 26; 11 27; 11 29; 11 39; 11 45; 12 7; 12 8; 12 13; 12 15; 12 16; 12 30; 12 31; 12 37; 12 43; 13 3; 13 4; 13 22; 13 37; 13 40; 13 42; 14 23; 14 25; 14 31; 14 46; 14 48; 15 19; 15 50; 15 51; 15 60; 16 2; 16 12; 16 23; 16 25; 16 31; 16 38; 17 2; 17 18; 17 33; 17 39; 17 54; 18 14; 18 15; 18 23; 18 25; 18 46; 19 20; 19 27; 19 30; 19 42; 19 55; 19 56; 19 58; 20 2; 20 8; 20 17; 20 29; 20 31; 20 40; 20 45; 20 49; 21 44; 21 45; 22 7; 22 15; 22 17; 22 18; 22 33; 22 45; 22 52; 22 54; 23 20; 23 50; 24 31; 24 51; 25 5; 25 7; 25 17; 25 31; 25 36; 25 41; 25 59; 26 6; 26 8; 26 11; 26 20; 26 21; 26 22; 26 35; 26 45; 26 48; 26 49; 26 57; 27 23; 27 45; 27 51; 28 13; 28 29; 28 51; 28 55; 28 59; 29 7; 29 23; 29 25; 30 3; 30 13; 30 15; 30 23; 30 49; 30 52; 31 43; 31 45; 31 53; 31 60; 32 3; 32 11; 32 30; 32 48; 32 49; 32 52; 32 55; 33 14; 33 16; 33 18; 33 19; 33 34; 33 39; 33 42; 33 44; 33 45; 33 59; 34 17; 34 36; 34 40; 34 48; 35 6; 35 13; 35 22; 35 24; 35 47; 35 55; 36 23; 36 27; 36 30; 36 33; 36 43; 36 50; 36 52; 37 14; 37 43; 37 53; 37 60; 38 3; 38 13; 38 57; 38 58; 39 15; 39 30; 39 34; 39 35; 39 40; 39 46; 39 49; 39 54; 40 8; 40 21; 40 29; 40 59; 41 10; 41 16; 41 39; 41 49; 41 56; 42 5; 42 28; 42 52; 42 57; 43 7; 43 19; 43 24; 44 3; 44 9; 44 14; 44 25; 44 28; 44 29; 45 8; 45 35; 45 39; 45 49; 45 50; 45 51; 45 57; 46 16; 46 36; 47 2; 47 15; 47 48; 47 55; 47 59; 48 9; 48 25; 48 29; 48 42; 48 52; 48 56; 49 8; 49 12; 49 33; 50 13; 50 19; 50 30; 50 52; 50 55; 51 10; 51 18; 51 28; 51 29; 51 36; 51 38; 51 43; 51 47; 52 2; 52 3; 52 4; 52 8; 52 23; 52 34; 53 4; 53 9; 53 14; 53 23; 53 33; 53 34; 53 35; 53 39; 53 51; 53 56; 53 59; 53 60; 54 4; 54 8; 54 45; 54 53; 55 6; 55 11; 55 19; 55 26; 55 29; 55 39; 56 14; 56 37; 56 39; 56 57; 57 21; 57 25; 57 26; 57 30; 57 40; 57 50; 58 4; 58 11; 58 16; 58 18; 58 21; 58 26; 58 45; 58 54; 59 12; 59 27; 59 28; 59 42; 59 45; 59 51]
global d_x = [1.0, 4.0, 1.0, 5.0, 10.0, 3.0, 7.0, 4.0, 6.0, 1.0, 7.0, 5.0, 8.0, 5.0, 10.0, 10.0, 7.0, 8.0, 6.0, 5.0, 4.0, 9.0, 1.0, 7.0, 1.0, 5.0, 8.0, 4.0, 2.0, 10.0, 9.0, 10.0, 10.0, 4.0, 6.0, 2.0, 2.0, 9.0, 9.0, 5.0, 2.0, 3.0, 5.0, 9.0, 4.0, 9.0, 3.0, 2.0, 6.0, 5.0, 4.0, 4.0, 10.0, 1.0, 9.0, 7.0, 1.0, 5.0, 6.0, 1.0, 4.0, 6.0, 2.0, 2.0, 2.0, 7.0, 4.0, 6.0, 5.0, 3.0, 2.0, 9.0, 1.0, 7.0, 4.0, 6.0, 3.0, 6.0, 7.0, 10.0, 2.0, 4.0, 10.0, 7.0, 9.0, 6.0, 10.0, 8.0, 6.0, 2.0, 7.0, 10.0, 9.0, 5.0, 7.0, 1.0, 9.0, 1.0, 7.0, 1.0, 6.0, 2.0, 6.0, 4.0, 6.0, 1.0, 2.0, 5.0, 3.0, 2.0, 1.0, 4.0, 2.0, 4.0, 7.0, 2.0, 6.0, 2.0, 2.0, 2.0, 6.0, 8.0, 2.0, 3.0, 1.0, 5.0, 5.0, 8.0, 4.0, 4.0, 3.0, 2.0, 5.0, 10.0, 3.0, 4.0, 1.0, 7.0, 1.0, 10.0, 8.0, 8.0, 5.0, 7.0, 6.0, 5.0, 5.0, 2.0, 9.0, 1.0, 2.0, 9.0, 2.0, 8.0, 4.0, 10.0, 5.0, 7.0, 8.0, 3.0, 10.0, 2.0, 8.0, 10.0, 9.0, 3.0, 8.0, 2.0, 5.0, 10.0, 6.0, 4.0, 1.0, 5.0, 10.0, 4.0, 6.0, 2.0, 1.0, 3.0, 1.0, 8.0, 2.0, 1.0, 7.0, 3.0, 5.0, 4.0, 3.0, 8.0, 10.0, 9.0, 3.0, 10.0, 1.0, 9.0, 10.0, 1.0, 3.0, 7.0, 5.0, 8.0, 7.0, 1.0, 5.0, 7.0, 8.0, 4.0, 1.0, 9.0, 5.0, 6.0, 5.0, 6.0, 5.0, 9.0, 3.0, 9.0, 7.0, 8.0, 4.0, 1.0, 3.0, 2.0, 1.0, 5.0, 1.0, 9.0, 8.0, 6.0, 5.0, 4.0, 9.0, 2.0, 4.0, 4.0, 5.0, 6.0, 10.0, 1.0, 7.0, 5.0, 6.0, 8.0, 10.0, 9.0, 1.0, 4.0, 10.0, 7.0, 4.0, 10.0, 1.0, 10.0, 1.0, 5.0, 10.0, 3.0, 9.0, 7.0, 6.0, 2.0, 2.0, 5.0, 3.0, 6.0, 5.0, 7.0, 5.0, 8.0, 3.0, 6.0, 1.0, 3.0, 8.0, 3.0, 10.0, 9.0, 7.0, 7.0, 4.0, 6.0, 4.0, 1.0, 1.0, 1.0, 9.0, 9.0, 1.0, 7.0, 3.0, 6.0, 10.0, 9.0, 9.0, 7.0, 4.0, 10.0, 3.0, 4.0, 10.0, 5.0, 9.0, 5.0, 8.0, 4.0, 3.0, 8.0, 8.0, 6.0, 5.0, 1.0, 9.0, 10.0, 6.0, 10.0, 1.0, 10.0, 3.0, 2.0, 6.0, 10.0, 4.0, 3.0, 3.0, 7.0, 6.0, 8.0, 10.0]
global b_x = 5
global d_y = [3.0, 6.0, 5.0, 2.0, 3.0, 10.0, 1.0, 9.0, 9.0, 9.0, 9.0, 9.0, 8.0, 6.0, 6.0, 8.0, 5.0, 2.0, 3.0, 3.0, 4.0, 1.0, 3.0, 6.0, 3.0, 7.0, 5.0, 4.0, 10.0, 4.0, 6.0, 6.0, 5.0, 7.0, 7.0, 9.0, 5.0, 2.0, 10.0, 8.0, 8.0, 3.0, 10.0, 10.0, 1.0, 4.0, 3.0, 1.0, 1.0, 1.0, 8.0, 2.0, 2.0, 1.0, 4.0, 2.0, 1.0, 5.0, 3.0, 10.0, 7.0, 9.0, 1.0, 10.0, 3.0, 6.0, 5.0, 10.0, 8.0, 3.0, 8.0, 6.0, 3.0, 3.0, 3.0, 7.0, 2.0, 6.0, 8.0, 5.0, 2.0, 5.0, 5.0, 5.0, 2.0, 5.0, 9.0, 8.0, 7.0, 7.0, 7.0, 2.0, 7.0, 1.0, 8.0, 8.0, 10.0, 7.0, 9.0, 8.0, 9.0, 6.0, 7.0, 5.0, 8.0, 9.0, 10.0, 1.0, 6.0, 8.0, 9.0, 10.0, 4.0, 6.0, 4.0, 6.0, 6.0, 9.0, 9.0, 7.0, 2.0, 1.0, 6.0, 3.0, 5.0, 1.0, 8.0, 5.0, 6.0, 7.0, 2.0, 7.0, 6.0, 8.0, 6.0, 9.0, 4.0, 5.0, 1.0, 1.0, 2.0, 4.0, 7.0, 10.0, 6.0, 9.0, 2.0, 4.0, 3.0, 6.0, 6.0, 9.0, 8.0, 3.0, 4.0, 5.0, 9.0, 6.0, 7.0, 10.0, 7.0, 3.0, 6.0, 8.0, 9.0, 9.0, 1.0, 5.0, 3.0, 9.0, 3.0, 8.0, 6.0, 7.0, 3.0, 8.0, 10.0, 9.0, 9.0, 5.0, 3.0, 5.0, 4.0, 7.0, 1.0, 1.0, 6.0, 7.0, 6.0, 6.0, 2.0, 9.0, 9.0, 4.0, 5.0, 8.0, 5.0, 7.0, 6.0, 1.0, 1.0, 5.0, 3.0, 3.0, 7.0, 2.0, 3.0, 2.0, 4.0, 6.0, 1.0, 2.0, 10.0, 7.0, 9.0, 6.0, 7.0, 8.0, 1.0, 7.0, 1.0, 4.0, 7.0, 2.0, 5.0, 6.0, 4.0, 5.0, 5.0, 5.0, 5.0, 2.0, 5.0, 5.0, 1.0, 1.0, 7.0, 3.0, 5.0, 2.0, 6.0, 4.0, 4.0, 10.0, 1.0, 1.0, 4.0, 6.0, 1.0, 1.0, 5.0, 7.0, 10.0, 9.0, 8.0, 10.0, 1.0, 5.0, 1.0, 10.0, 7.0, 4.0, 10.0, 7.0, 8.0, 8.0, 10.0, 3.0, 1.0, 5.0, 10.0, 5.0, 3.0, 3.0, 6.0, 10.0, 5.0, 4.0, 6.0, 4.0, 7.0, 4.0, 10.0, 9.0, 8.0, 6.0, 3.0, 3.0, 8.0, 9.0, 7.0, 4.0, 9.0, 4.0, 2.0, 4.0, 4.0, 1.0, 9.0, 2.0, 8.0, 1.0, 9.0, 7.0, 9.0, 8.0, 4.0, 10.0, 5.0, 10.0, 10.0, 10.0, 2.0, 2.0, 8.0, 10.0, 6.0, 6.0, 3.0, 7.0, 1.0, 5.0, 2.0, 3.0, 1.0, 3.0, 3.0, 9.0, 2.0]
global b_y = 10
global p = [0.839, 0.182, 0.353, 0.527, 0.142, 0.343, 0.813, 0.967, 0.42, 0.78, 0.5, 0.974, 0.864, 0.105, 0.282, 0.814, 0.823, 0.559, 0.519, 0.428, 0.085, 0.45, 0.329, 0.173, 0.356, 0.247, 0.118, 0.092, 0.497, 0.588, 0.899, 0.48, 0.285, 0.682, 0.668, 0.351, 0.659, 0.223, 0.37, 0.715, 0.783, 0.647, 0.924, 0.545, 0.982, 0.373, 0.932, 0.636, 0.025, 0.296, 0.6, 0.975, 0.499, 0.123, 0.61, 0.513, 0.998, 0.317, 0.07, 0.191, 0.113, 0.691, 0.333, 0.806, 0.074, 0.633, 0.343, 0.62, 0.02, 0.762, 0.239, 0.702, 0.381, 0.79, 0.622, 0.468, 0.624, 0.466, 0.666, 0.662, 0.421, 0.189, 0.94, 0.308, 0.269, 0.4, 0.129, 0.623, 0.661, 0.511, 0.479, 0.336, 0.013, 0.415, 0.852, 0.779, 0.03, 0.357, 0.512, 0.201, 0.169, 0.85, 0.921, 0.025, 0.797, 0.523, 0.325, 0.39, 0.473, 0.529, 0.066, 0.359, 0.909, 0.738, 0.949, 0.785, 0.771, 0.868, 0.666, 0.096, 0.185, 0.362, 0.678, 0.075, 0.451, 0.402, 0.134, 0.407, 0.21, 0.595, 0.085, 0.48, 0.663, 0.282, 0.151, 0.401, 0.168, 0.442, 0.072, 0.023, 0.275, 0.029, 0.149, 0.588, 0.61, 0.62, 0.144, 0.291, 0.526, 0.676, 0.609, 0.476, 0.032, 0.45, 0.752, 0.698, 0.183, 0.695, 0.888, 0.862, 0.293, 0.028, 0.193, 0.986, 0.202, 0.781, 0.902, 0.418, 0.458, 0.955, 0.176, 0.6, 0.59, 0.286, 0.973, 0.169, 0.239, 0.857, 0.213, 0.541, 0.607, 0.059, 0.333, 0.403, 0.784, 0.567, 0.161, 0.619, 0.203, 0.263, 0.187, 0.568, 0.596, 0.048, 0.116, 0.679, 0.581, 0.96, 0.941, 0.432, 0.912, 0.299, 0.967, 0.724, 0.31, 0.824, 0.765, 0.89, 0.209, 0.921, 0.057, 0.972, 0.355, 0.611, 0.352, 0.037, 0.499, 0.029, 0.773, 0.463, 0.014, 0.744, 0.078, 0.765, 0.632, 0.891, 0.107, 0.728, 0.703, 0.316, 0.553, 0.738, 0.009, 0.042, 0.775, 0.563, 0.209, 0.628, 0.161, 0.965, 0.509, 0.508, 0.985, 0.208, 0.784, 0.805, 0.86, 0.309, 0.337, 0.633, 0.588, 0.891, 0.752, 0.03, 0.128, 0.536, 0.574, 0.936, 0.224, 0.438, 0.258, 0.947, 0.163, 0.553, 0.707, 0.789, 0.437, 0.249, 0.842, 0.29, 0.042, 0.22, 0.53, 0.694, 0.245, 0.878, 0.127, 0.67, 0.285, 0.849, 0.934, 0.846, 0.503, 0.858, 0.438, 0.683, 0.568, 0.647, 0.677, 0.242, 0.105, 0.352, 0.014, 0.902, 0.517, 0.736, 0.981, 0.976, 0.946, 0.728, 0.293, 0.297, 0.514, 0.601, 0.511, 0.607, 0.649, 0.9, 0.326, 0.214, 0.597, 0.499, 0.5, 0.574, 0.277, 0.643, 0.212, 0.269, 0.727, 0.27, 0.012, 0.014, 0.763, 0.828, 0.016, 0.87, 0.537, 0.483, 0.942]
global q = [0.864, 0.209, 0.69, 0.981, 0.255, 0.68, 0.844, 0.967, 0.85, 0.889, 0.708, 0.983, 0.883, 0.191, 0.485, 0.932, 0.979, 0.625, 0.624, 0.597, 0.516, 0.774, 0.373, 0.665, 0.45, 0.634, 0.214, 0.582, 0.508, 0.745, 0.941, 0.605, 0.381, 0.86, 0.949, 0.919, 0.837, 0.839, 0.61, 0.726, 0.823, 0.704, 0.998, 0.553, 0.982, 0.928, 0.953, 0.926, 0.902, 0.913, 0.952, 0.99, 0.865, 0.127, 0.926, 0.783, 0.999, 0.766, 0.637, 0.894, 0.655, 0.841, 0.942, 0.985, 0.669, 0.847, 0.358, 0.87, 0.211, 0.79, 0.298, 0.906, 0.594, 0.951, 0.644, 0.606, 0.682, 0.753, 0.798, 0.859, 0.987, 0.256, 0.994, 0.723, 0.7, 0.633, 0.918, 0.731, 0.817, 0.658, 0.857, 0.608, 0.475, 0.822, 0.986, 0.898, 0.458, 0.946, 0.598, 0.253, 0.292, 0.969, 0.959, 0.368, 0.972, 0.759, 0.493, 0.825, 0.788, 0.705, 0.622, 0.518, 0.984, 0.75, 0.998, 0.984, 0.917, 0.945, 0.972, 0.358, 0.347, 0.678, 0.817, 0.184, 0.686, 0.965, 0.706, 0.817, 0.939, 0.877, 0.904, 0.99, 0.978, 0.831, 0.402, 0.425, 0.554, 0.567, 0.263, 0.189, 0.627, 0.567, 0.826, 0.931, 0.644, 0.696, 0.427, 0.611, 0.552, 0.932, 0.828, 0.926, 0.498, 0.851, 0.928, 0.934, 0.683, 0.869, 0.995, 0.988, 0.335, 0.63, 0.959, 0.989, 0.669, 0.895, 0.966, 0.545, 0.61, 0.959, 0.296, 0.832, 0.921, 0.811, 0.981, 0.366, 0.752, 0.86, 0.931, 0.555, 0.843, 0.513, 0.743, 0.902, 0.888, 0.931, 0.646, 0.998, 0.419, 0.911, 0.217, 0.696, 0.778, 0.368, 0.94, 0.929, 0.608, 0.978, 0.984, 0.681, 0.937, 0.824, 0.97, 0.73, 0.597, 0.995, 0.992, 0.898, 0.854, 0.997, 0.445, 0.997, 0.844, 0.67, 0.908, 0.404, 0.62, 0.798, 0.951, 0.607, 0.046, 0.862, 0.404, 0.828, 0.936, 0.937, 0.697, 0.811, 0.716, 0.486, 0.958, 0.845, 0.145, 0.043, 0.781, 0.764, 0.301, 0.811, 0.23, 0.998, 0.846, 0.943, 0.988, 0.353, 0.93, 0.992, 0.999, 0.874, 0.915, 0.757, 0.683, 0.99, 0.819, 0.306, 0.417, 0.542, 0.835, 0.947, 0.832, 0.537, 0.69, 0.965, 0.775, 0.742, 0.898, 0.86, 0.77, 0.46, 0.87, 0.794, 0.384, 0.594, 0.69, 0.99, 0.739, 0.956, 0.133, 0.693, 0.735, 0.85, 0.993, 0.969, 0.686, 0.967, 0.817, 0.744, 0.737, 0.843, 0.843, 0.391, 0.59, 0.813, 0.771, 0.978, 0.688, 0.874, 0.998, 0.983, 0.955, 0.796, 0.454, 0.753, 0.959, 0.655, 0.723, 0.705, 0.717, 0.934, 0.38, 0.672, 0.768, 0.812, 0.746, 0.759, 0.526, 0.992, 0.783, 0.894, 0.813, 0.811, 0.568, 0.148, 0.785, 0.999, 0.609, 0.917, 0.888, 0.862, 0.959]
global origin = 1
global destination = 60