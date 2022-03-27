global arcs = [1 14; 1 16; 1 28; 2 12; 2 15; 2 25; 2 34; 2 40; 2 58; 3 34; 3 58; 4 2; 4 8; 4 18; 4 21; 4 32; 4 42; 5 3; 5 11; 5 23; 5 37; 5 44; 5 46; 5 48; 6 15; 6 18; 6 19; 6 20; 6 32; 6 40; 6 59; 7 8; 7 13; 7 21; 7 34; 7 55; 7 57; 8 5; 8 9; 8 10; 8 21; 8 26; 8 37; 8 42; 8 56; 9 7; 9 10; 9 26; 9 46; 9 57; 10 3; 10 4; 10 14; 10 17; 10 18; 10 19; 10 22; 10 31; 10 40; 10 46; 10 55; 10 60; 11 28; 11 36; 11 40; 11 56; 12 2; 12 15; 12 20; 12 28; 12 32; 12 49; 12 55; 12 56; 13 3; 13 6; 13 7; 13 51; 13 54; 13 58; 14 6; 14 10; 14 17; 14 31; 14 40; 14 43; 14 52; 15 12; 15 30; 15 51; 15 54; 16 13; 16 23; 16 43; 16 45; 16 50; 16 58; 17 7; 17 9; 17 36; 17 46; 17 54; 18 7; 18 8; 18 31; 18 41; 19 12; 19 45; 19 60; 20 4; 20 26; 20 38; 20 41; 20 51; 21 5; 21 46; 22 4; 22 14; 22 24; 22 29; 22 41; 22 49; 22 52; 23 15; 23 28; 23 35; 23 43; 23 51; 23 57; 24 12; 24 21; 24 29; 24 42; 25 8; 25 19; 25 32; 26 18; 26 21; 26 37; 26 43; 26 46; 26 59; 27 10; 27 26; 27 28; 27 29; 27 33; 28 5; 28 25; 28 48; 28 58; 29 10; 29 15; 29 19; 29 21; 29 41; 29 44; 29 51; 29 53; 29 54; 30 20; 30 24; 30 34; 30 39; 30 43; 31 5; 31 16; 31 17; 31 22; 31 34; 31 55; 32 22; 32 33; 33 18; 33 22; 33 43; 33 51; 33 54; 34 12; 34 16; 34 23; 34 38; 34 39; 34 49; 35 3; 35 9; 35 29; 35 53; 36 16; 36 21; 36 31; 36 40; 36 48; 36 53; 36 54; 37 14; 37 22; 37 24; 37 40; 37 41; 37 43; 37 52; 37 58; 38 15; 38 23; 38 28; 38 33; 38 51; 38 52; 39 18; 39 24; 39 56; 40 6; 40 9; 40 10; 40 24; 40 51; 40 53; 41 14; 41 27; 41 37; 41 38; 41 40; 42 3; 42 6; 42 19; 42 22; 42 26; 42 39; 42 45; 42 50; 42 59; 42 60; 43 6; 43 11; 43 18; 43 21; 44 3; 44 12; 44 13; 44 30; 44 51; 44 52; 45 3; 45 22; 45 57; 45 60; 46 13; 46 25; 46 31; 46 39; 46 58; 46 60; 47 19; 47 58; 48 3; 48 14; 48 28; 48 29; 48 37; 48 40; 48 44; 48 59; 49 10; 49 26; 49 29; 49 39; 49 47; 49 51; 49 53; 49 55; 50 35; 50 40; 50 46; 50 55; 51 5; 51 10; 51 12; 51 33; 51 38; 51 41; 51 43; 52 10; 52 15; 52 27; 52 39; 52 44; 52 48; 53 4; 53 10; 53 11; 53 42; 53 49; 53 54; 54 6; 54 8; 54 23; 54 24; 54 29; 54 45; 54 58; 55 4; 55 9; 55 35; 55 46; 55 48; 55 53; 55 57; 56 30; 56 36; 56 45; 56 52; 57 14; 57 32; 57 41; 57 51; 57 58; 58 9; 58 10; 58 31; 58 34; 59 13; 59 25; 59 28; 59 45; 59 47; 59 53; 59 54]
global d_x = [6.0, 2.0, 9.0, 1.0, 1.0, 6.0, 2.0, 9.0, 9.0, 7.0, 5.0, 1.0, 3.0, 5.0, 4.0, 3.0, 1.0, 2.0, 7.0, 6.0, 3.0, 10.0, 2.0, 10.0, 2.0, 10.0, 7.0, 10.0, 2.0, 8.0, 6.0, 6.0, 8.0, 8.0, 10.0, 2.0, 2.0, 9.0, 5.0, 1.0, 6.0, 1.0, 1.0, 9.0, 9.0, 5.0, 6.0, 6.0, 6.0, 5.0, 5.0, 4.0, 7.0, 5.0, 6.0, 1.0, 10.0, 6.0, 1.0, 4.0, 3.0, 10.0, 4.0, 2.0, 6.0, 9.0, 9.0, 8.0, 4.0, 4.0, 8.0, 6.0, 1.0, 1.0, 6.0, 3.0, 6.0, 4.0, 9.0, 10.0, 3.0, 2.0, 2.0, 4.0, 10.0, 3.0, 1.0, 7.0, 9.0, 10.0, 7.0, 6.0, 6.0, 2.0, 1.0, 4.0, 8.0, 3.0, 10.0, 8.0, 3.0, 5.0, 7.0, 7.0, 5.0, 3.0, 9.0, 9.0, 9.0, 5.0, 9.0, 2.0, 7.0, 8.0, 8.0, 8.0, 4.0, 3.0, 3.0, 4.0, 5.0, 10.0, 6.0, 1.0, 4.0, 10.0, 3.0, 10.0, 2.0, 1.0, 6.0, 7.0, 7.0, 6.0, 2.0, 5.0, 1.0, 2.0, 4.0, 8.0, 8.0, 6.0, 6.0, 8.0, 4.0, 5.0, 8.0, 6.0, 3.0, 4.0, 5.0, 8.0, 5.0, 5.0, 7.0, 10.0, 5.0, 7.0, 9.0, 8.0, 2.0, 5.0, 4.0, 5.0, 10.0, 5.0, 10.0, 8.0, 8.0, 5.0, 8.0, 1.0, 10.0, 6.0, 10.0, 4.0, 7.0, 1.0, 6.0, 3.0, 2.0, 3.0, 4.0, 9.0, 3.0, 10.0, 8.0, 4.0, 4.0, 10.0, 10.0, 3.0, 9.0, 9.0, 6.0, 8.0, 1.0, 7.0, 3.0, 7.0, 5.0, 1.0, 2.0, 7.0, 4.0, 6.0, 9.0, 5.0, 2.0, 9.0, 7.0, 6.0, 10.0, 10.0, 6.0, 6.0, 3.0, 2.0, 7.0, 6.0, 10.0, 4.0, 9.0, 2.0, 10.0, 5.0, 1.0, 9.0, 9.0, 2.0, 5.0, 6.0, 8.0, 7.0, 1.0, 9.0, 5.0, 6.0, 10.0, 3.0, 10.0, 8.0, 8.0, 6.0, 10.0, 10.0, 7.0, 7.0, 4.0, 9.0, 10.0, 5.0, 2.0, 7.0, 5.0, 3.0, 1.0, 10.0, 9.0, 9.0, 6.0, 5.0, 7.0, 9.0, 9.0, 7.0, 3.0, 9.0, 10.0, 10.0, 9.0, 9.0, 8.0, 10.0, 5.0, 8.0, 2.0, 7.0, 3.0, 2.0, 4.0, 1.0, 7.0, 6.0, 2.0, 2.0, 9.0, 6.0, 5.0, 4.0, 10.0, 7.0, 9.0, 4.0, 10.0, 1.0, 9.0, 2.0, 2.0, 2.0, 8.0, 9.0, 1.0, 5.0, 7.0, 2.0, 2.0, 7.0, 3.0, 5.0, 10.0, 1.0, 5.0, 7.0, 4.0, 3.0, 5.0, 9.0, 9.0, 2.0, 2.0, 7.0, 7.0, 6.0, 5.0, 3.0, 10.0, 8.0]
global b_x = 5
global d_y = [3.0, 6.0, 7.0, 2.0, 8.0, 6.0, 7.0, 1.0, 3.0, 3.0, 7.0, 9.0, 1.0, 3.0, 2.0, 8.0, 9.0, 1.0, 3.0, 7.0, 8.0, 5.0, 5.0, 3.0, 6.0, 10.0, 9.0, 2.0, 4.0, 3.0, 8.0, 10.0, 1.0, 1.0, 8.0, 3.0, 6.0, 7.0, 1.0, 9.0, 8.0, 1.0, 7.0, 6.0, 1.0, 5.0, 4.0, 2.0, 3.0, 9.0, 6.0, 2.0, 1.0, 1.0, 5.0, 3.0, 3.0, 9.0, 7.0, 1.0, 2.0, 7.0, 9.0, 1.0, 6.0, 2.0, 9.0, 3.0, 7.0, 6.0, 7.0, 4.0, 2.0, 9.0, 1.0, 7.0, 6.0, 8.0, 2.0, 6.0, 1.0, 9.0, 10.0, 7.0, 4.0, 5.0, 8.0, 7.0, 1.0, 1.0, 4.0, 6.0, 5.0, 9.0, 2.0, 9.0, 1.0, 3.0, 2.0, 3.0, 8.0, 10.0, 8.0, 8.0, 9.0, 4.0, 2.0, 7.0, 7.0, 9.0, 9.0, 6.0, 7.0, 10.0, 1.0, 8.0, 1.0, 9.0, 10.0, 5.0, 5.0, 7.0, 4.0, 4.0, 1.0, 6.0, 9.0, 3.0, 7.0, 4.0, 6.0, 10.0, 6.0, 2.0, 1.0, 2.0, 3.0, 6.0, 3.0, 9.0, 3.0, 2.0, 8.0, 9.0, 8.0, 3.0, 9.0, 4.0, 8.0, 5.0, 8.0, 9.0, 10.0, 8.0, 6.0, 1.0, 6.0, 8.0, 1.0, 8.0, 9.0, 4.0, 2.0, 2.0, 5.0, 8.0, 2.0, 2.0, 10.0, 3.0, 4.0, 9.0, 3.0, 9.0, 5.0, 3.0, 7.0, 2.0, 6.0, 8.0, 1.0, 3.0, 5.0, 10.0, 10.0, 10.0, 3.0, 4.0, 7.0, 2.0, 3.0, 6.0, 10.0, 4.0, 7.0, 7.0, 1.0, 6.0, 2.0, 4.0, 8.0, 4.0, 8.0, 2.0, 2.0, 3.0, 7.0, 4.0, 10.0, 10.0, 5.0, 1.0, 4.0, 5.0, 10.0, 5.0, 2.0, 9.0, 6.0, 4.0, 7.0, 2.0, 7.0, 7.0, 10.0, 5.0, 10.0, 3.0, 1.0, 6.0, 1.0, 3.0, 8.0, 7.0, 2.0, 5.0, 1.0, 8.0, 4.0, 8.0, 5.0, 1.0, 9.0, 10.0, 2.0, 7.0, 2.0, 8.0, 3.0, 7.0, 2.0, 1.0, 5.0, 8.0, 3.0, 7.0, 1.0, 7.0, 10.0, 4.0, 8.0, 7.0, 5.0, 3.0, 2.0, 8.0, 7.0, 6.0, 2.0, 7.0, 4.0, 4.0, 1.0, 6.0, 4.0, 8.0, 3.0, 2.0, 5.0, 5.0, 9.0, 3.0, 1.0, 6.0, 4.0, 1.0, 4.0, 8.0, 6.0, 6.0, 8.0, 8.0, 10.0, 1.0, 7.0, 6.0, 8.0, 3.0, 5.0, 1.0, 10.0, 6.0, 4.0, 6.0, 2.0, 10.0, 1.0, 9.0, 8.0, 9.0, 1.0, 1.0, 6.0, 7.0, 3.0, 7.0, 1.0, 8.0, 10.0, 9.0, 1.0, 4.0, 2.0, 8.0, 10.0, 9.0, 5.0, 2.0]
global b_y = 10
global p = [0.138, 0.672, 0.392, 0.82, 0.987, 0.184, 0.453, 0.99, 0.66, 0.139, 0.048, 0.637, 0.464, 0.68, 0.005, 0.405, 0.575, 0.788, 0.842, 0.278, 0.881, 0.378, 0.204, 0.191, 0.918, 0.543, 0.604, 0.357, 0.489, 0.686, 0.17, 0.171, 0.949, 0.188, 0.566, 0.463, 0.556, 0.734, 0.268, 0.984, 0.589, 0.787, 0.15, 0.716, 0.063, 0.589, 0.701, 0.683, 0.27, 0.595, 0.468, 0.652, 0.621, 0.212, 0.328, 0.164, 0.798, 0.867, 0.866, 0.697, 0.065, 0.247, 0.484, 0.25, 0.643, 0.985, 0.315, 0.966, 0.934, 0.524, 0.739, 0.093, 0.195, 0.64, 0.844, 0.721, 0.589, 0.028, 0.643, 0.882, 0.514, 0.988, 0.955, 0.121, 0.401, 0.193, 0.151, 0.908, 0.713, 0.577, 0.948, 0.209, 0.818, 0.223, 0.925, 0.07, 0.005, 0.985, 0.506, 0.775, 0.047, 0.663, 0.84, 0.615, 0.284, 0.856, 0.621, 0.514, 0.786, 0.249, 0.67, 0.945, 0.831, 0.735, 0.49, 0.326, 0.736, 0.715, 0.813, 0.465, 0.631, 0.704, 0.79, 0.625, 0.742, 0.936, 0.704, 0.57, 0.872, 0.41, 0.301, 0.784, 0.292, 0.445, 0.253, 0.951, 0.421, 0.424, 0.675, 0.216, 0.736, 0.61, 0.499, 0.113, 0.668, 0.381, 0.393, 0.736, 0.1, 0.368, 0.149, 0.855, 0.628, 0.304, 0.792, 0.686, 0.032, 0.075, 0.393, 0.695, 0.139, 0.19, 0.636, 0.184, 0.21, 0.695, 0.242, 0.391, 0.523, 0.792, 0.814, 0.178, 0.693, 0.688, 0.073, 0.433, 0.927, 0.391, 0.897, 0.965, 0.044, 0.768, 0.874, 0.789, 0.281, 0.773, 0.8, 0.738, 0.485, 0.321, 0.45, 0.278, 0.129, 0.464, 0.942, 0.168, 0.723, 0.263, 0.465, 0.609, 0.391, 0.006, 0.728, 0.318, 0.877, 0.547, 0.779, 0.75, 0.732, 0.71, 0.774, 0.675, 0.608, 0.053, 0.784, 0.714, 0.041, 0.452, 0.668, 0.901, 0.058, 0.214, 0.814, 0.414, 0.054, 0.444, 0.232, 0.413, 0.186, 0.104, 0.652, 0.321, 0.074, 0.51, 0.625, 0.512, 0.137, 0.654, 0.047, 0.549, 0.277, 0.096, 0.972, 0.788, 0.002, 0.477, 0.722, 0.327, 0.973, 0.307, 0.029, 0.745, 0.002, 0.055, 0.934, 0.516, 0.956, 0.484, 0.176, 0.975, 0.678, 0.097, 0.676, 0.243, 0.495, 0.92, 0.012, 0.403, 0.684, 0.561, 0.405, 0.712, 0.496, 0.035, 0.069, 0.191, 0.126, 0.468, 0.193, 0.871, 0.788, 0.989, 0.25, 0.274, 0.125, 0.559, 0.797, 0.48, 0.694, 0.421, 0.318, 0.929, 0.686, 0.547, 0.521, 0.703, 0.215, 0.524, 0.762, 0.037, 0.962, 0.37, 0.466, 0.133, 0.912, 0.402, 0.771, 0.261, 0.837, 0.516, 0.813, 0.594, 0.547, 0.274, 0.454, 0.308, 0.798, 0.918, 0.025, 0.872, 0.465, 0.891, 0.361, 0.042, 0.984, 0.799, 0.488, 0.372]
global q = [0.56, 0.911, 0.973, 0.996, 0.994, 0.285, 0.731, 0.992, 0.997, 0.891, 0.348, 0.945, 0.983, 0.906, 0.108, 0.42, 0.849, 0.888, 0.961, 0.657, 0.978, 0.555, 0.873, 0.196, 0.976, 0.798, 0.997, 0.958, 0.797, 0.951, 0.714, 0.489, 0.964, 0.651, 0.718, 0.632, 0.637, 0.99, 0.525, 0.986, 0.734, 0.876, 0.797, 0.806, 0.589, 0.715, 0.805, 0.888, 0.561, 0.598, 0.749, 0.854, 0.883, 0.637, 0.755, 0.633, 0.833, 0.896, 0.92, 0.982, 0.108, 0.708, 0.551, 0.926, 0.716, 0.996, 0.783, 0.986, 0.976, 0.804, 0.753, 0.97, 0.734, 0.93, 0.875, 0.894, 0.645, 0.891, 0.973, 0.916, 0.659, 0.991, 0.986, 0.571, 0.833, 0.773, 0.788, 0.911, 0.726, 0.949, 0.962, 0.827, 0.846, 0.299, 0.938, 0.668, 0.619, 0.985, 0.945, 0.81, 0.199, 0.754, 0.901, 0.877, 0.495, 0.858, 0.721, 0.683, 0.809, 0.903, 0.698, 0.977, 0.832, 0.889, 0.631, 0.954, 0.916, 0.789, 0.851, 0.826, 0.947, 0.878, 0.865, 0.638, 0.957, 0.991, 0.851, 0.767, 0.927, 0.413, 0.48, 0.973, 0.597, 0.596, 0.267, 0.97, 0.804, 0.948, 0.847, 0.78, 0.956, 0.676, 0.953, 0.36, 0.967, 0.467, 0.419, 0.95, 0.438, 0.564, 0.487, 0.895, 0.654, 0.928, 0.803, 0.836, 0.851, 0.799, 0.484, 0.912, 0.803, 0.517, 0.711, 0.657, 0.664, 0.767, 0.497, 0.76, 0.669, 0.947, 0.935, 0.562, 0.86, 0.857, 0.818, 0.946, 0.946, 0.438, 0.944, 0.973, 0.898, 0.873, 0.942, 0.892, 0.359, 0.99, 0.91, 0.851, 0.991, 0.653, 0.563, 0.387, 0.615, 0.897, 0.965, 0.342, 0.948, 0.386, 0.978, 0.864, 0.417, 0.732, 0.827, 0.91, 0.973, 0.906, 0.93, 0.892, 0.814, 0.943, 0.928, 0.902, 0.886, 0.992, 0.886, 0.731, 0.825, 0.454, 0.856, 0.918, 0.716, 0.946, 0.911, 0.517, 0.418, 0.668, 0.452, 0.566, 0.332, 0.515, 0.715, 0.919, 0.241, 0.75, 0.987, 0.833, 0.406, 0.97, 0.926, 0.597, 0.569, 0.166, 0.98, 0.912, 0.047, 0.655, 0.902, 0.668, 0.999, 0.946, 0.421, 0.98, 0.623, 0.653, 0.976, 0.717, 0.975, 0.719, 0.545, 0.984, 0.733, 0.97, 0.872, 0.586, 0.798, 0.968, 0.73, 0.705, 0.829, 0.759, 0.891, 0.82, 0.927, 0.986, 0.511, 0.294, 0.155, 0.693, 0.58, 0.979, 0.817, 0.999, 0.47, 0.694, 0.22, 0.986, 0.956, 0.823, 0.739, 0.64, 0.673, 0.989, 0.763, 0.701, 0.887, 0.825, 0.495, 0.657, 0.927, 0.203, 0.965, 0.676, 0.48, 0.364, 0.969, 0.724, 0.877, 0.34, 0.989, 0.798, 0.814, 0.697, 0.783, 0.419, 0.751, 0.486, 0.91, 0.987, 0.232, 0.926, 0.92, 0.973, 0.662, 0.326, 0.993, 0.941, 0.71, 0.973]
global origin = 1
global destination = 60