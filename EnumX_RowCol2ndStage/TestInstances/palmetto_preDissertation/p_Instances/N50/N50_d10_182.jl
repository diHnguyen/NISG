global arcs = [1 2; 1 15; 1 16; 1 17; 1 23; 1 24; 1 31; 1 42; 2 9; 2 12; 2 33; 2 43; 3 8; 3 24; 3 41; 3 48; 4 2; 4 14; 4 19; 4 21; 4 30; 4 32; 4 37; 5 3; 5 13; 5 36; 5 47; 5 50; 6 15; 6 16; 6 22; 6 26; 7 9; 7 21; 7 26; 7 30; 7 40; 7 42; 7 46; 8 3; 8 5; 8 30; 8 31; 9 3; 9 12; 9 15; 9 39; 10 12; 10 29; 10 33; 10 39; 10 47; 11 14; 11 30; 11 38; 12 11; 12 20; 12 28; 12 39; 13 4; 13 8; 13 16; 13 30; 13 37; 13 39; 14 16; 14 23; 14 26; 14 40; 14 50; 15 2; 15 9; 15 30; 15 33; 16 9; 16 10; 16 22; 16 23; 16 24; 16 29; 16 34; 16 41; 16 42; 16 43; 16 47; 17 4; 17 24; 17 30; 17 34; 17 47; 18 15; 18 24; 18 33; 18 36; 18 42; 18 44; 19 3; 19 17; 19 22; 19 30; 19 35; 19 40; 19 43; 20 8; 20 14; 20 33; 21 5; 21 13; 21 15; 21 25; 21 35; 22 12; 22 18; 22 21; 22 35; 22 49; 23 3; 23 20; 23 22; 23 24; 23 28; 23 31; 24 12; 24 19; 25 5; 25 24; 25 27; 25 38; 25 39; 26 4; 26 34; 26 37; 26 40; 26 41; 27 2; 27 36; 27 37; 28 5; 28 24; 28 26; 28 36; 28 46; 28 48; 29 4; 29 8; 29 16; 29 25; 29 34; 29 47; 30 4; 30 11; 30 19; 30 22; 30 27; 30 34; 30 46; 31 2; 31 22; 31 26; 31 29; 31 37; 32 6; 32 19; 32 50; 33 7; 33 12; 33 14; 33 16; 33 25; 33 27; 33 37; 33 38; 33 43; 33 45; 33 50; 34 13; 34 23; 34 30; 34 40; 34 43; 34 44; 34 45; 34 48; 35 5; 35 20; 35 29; 35 34; 35 43; 35 44; 36 13; 36 19; 36 30; 36 43; 37 2; 37 7; 37 16; 37 21; 37 30; 37 38; 37 39; 37 43; 38 6; 38 30; 38 34; 38 49; 39 2; 39 18; 39 26; 39 27; 39 34; 39 38; 39 45; 39 50; 40 5; 40 11; 40 17; 40 21; 40 37; 40 39; 40 42; 41 13; 41 14; 41 15; 41 16; 41 35; 41 48; 42 15; 42 21; 42 25; 42 41; 42 47; 43 28; 43 35; 43 44; 44 5; 44 15; 44 36; 44 42; 44 48; 45 7; 45 11; 45 14; 45 24; 45 30; 45 36; 45 47; 46 29; 46 38; 46 47; 47 4; 47 19; 47 32; 47 37; 47 49; 48 3; 48 17; 48 35; 49 8; 49 23; 49 29; 49 34; 49 44]
global d_x = [6.0, 3.0, 4.0, 7.0, 3.0, 5.0, 3.0, 2.0, 6.0, 4.0, 8.0, 2.0, 1.0, 9.0, 4.0, 7.0, 4.0, 9.0, 2.0, 6.0, 9.0, 3.0, 1.0, 2.0, 5.0, 1.0, 9.0, 1.0, 5.0, 7.0, 3.0, 7.0, 5.0, 3.0, 4.0, 1.0, 8.0, 7.0, 2.0, 3.0, 3.0, 6.0, 9.0, 9.0, 10.0, 8.0, 2.0, 4.0, 4.0, 6.0, 2.0, 5.0, 7.0, 7.0, 7.0, 2.0, 10.0, 9.0, 2.0, 6.0, 8.0, 3.0, 5.0, 3.0, 8.0, 2.0, 8.0, 5.0, 7.0, 8.0, 1.0, 9.0, 6.0, 8.0, 6.0, 1.0, 4.0, 1.0, 1.0, 1.0, 4.0, 2.0, 1.0, 10.0, 2.0, 5.0, 8.0, 6.0, 7.0, 8.0, 10.0, 2.0, 8.0, 9.0, 8.0, 1.0, 9.0, 10.0, 7.0, 5.0, 8.0, 9.0, 4.0, 3.0, 9.0, 7.0, 8.0, 7.0, 10.0, 1.0, 6.0, 6.0, 2.0, 2.0, 4.0, 6.0, 3.0, 9.0, 6.0, 7.0, 5.0, 4.0, 3.0, 10.0, 1.0, 4.0, 7.0, 5.0, 9.0, 6.0, 10.0, 4.0, 7.0, 8.0, 10.0, 6.0, 6.0, 8.0, 6.0, 10.0, 9.0, 9.0, 1.0, 4.0, 10.0, 9.0, 5.0, 6.0, 10.0, 8.0, 6.0, 6.0, 10.0, 9.0, 1.0, 2.0, 6.0, 2.0, 2.0, 6.0, 7.0, 1.0, 2.0, 9.0, 3.0, 8.0, 3.0, 3.0, 7.0, 8.0, 10.0, 9.0, 5.0, 3.0, 1.0, 10.0, 10.0, 4.0, 7.0, 7.0, 6.0, 1.0, 4.0, 4.0, 5.0, 2.0, 5.0, 7.0, 2.0, 10.0, 4.0, 3.0, 8.0, 7.0, 1.0, 2.0, 8.0, 2.0, 5.0, 1.0, 3.0, 5.0, 9.0, 4.0, 1.0, 10.0, 9.0, 6.0, 4.0, 5.0, 9.0, 10.0, 8.0, 3.0, 7.0, 7.0, 8.0, 5.0, 5.0, 5.0, 2.0, 6.0, 6.0, 7.0, 5.0, 9.0, 8.0, 4.0, 3.0, 7.0, 7.0, 2.0, 2.0, 9.0, 9.0, 10.0, 10.0, 8.0, 3.0, 7.0, 3.0, 4.0, 2.0, 8.0, 5.0, 4.0, 4.0, 5.0, 1.0, 3.0, 1.0, 5.0, 4.0, 1.0, 5.0, 7.0, 6.0, 6.0, 9.0, 7.0, 2.0, 2.0]
global b_x = 5
global d_y = [2.0, 8.0, 5.0, 10.0, 7.0, 10.0, 5.0, 5.0, 3.0, 7.0, 2.0, 10.0, 8.0, 5.0, 9.0, 8.0, 7.0, 3.0, 1.0, 2.0, 6.0, 6.0, 6.0, 2.0, 3.0, 6.0, 10.0, 9.0, 10.0, 1.0, 8.0, 6.0, 3.0, 8.0, 5.0, 6.0, 6.0, 6.0, 4.0, 6.0, 7.0, 4.0, 7.0, 2.0, 6.0, 7.0, 4.0, 8.0, 4.0, 1.0, 1.0, 1.0, 10.0, 3.0, 2.0, 2.0, 7.0, 7.0, 8.0, 5.0, 5.0, 10.0, 10.0, 6.0, 5.0, 4.0, 4.0, 4.0, 4.0, 4.0, 8.0, 10.0, 4.0, 10.0, 3.0, 3.0, 8.0, 2.0, 2.0, 1.0, 6.0, 7.0, 3.0, 3.0, 6.0, 8.0, 6.0, 2.0, 2.0, 5.0, 3.0, 5.0, 10.0, 10.0, 6.0, 10.0, 8.0, 1.0, 6.0, 2.0, 9.0, 7.0, 9.0, 10.0, 9.0, 4.0, 10.0, 4.0, 6.0, 1.0, 6.0, 7.0, 9.0, 4.0, 3.0, 9.0, 8.0, 7.0, 9.0, 8.0, 7.0, 7.0, 1.0, 1.0, 5.0, 9.0, 3.0, 5.0, 10.0, 6.0, 4.0, 4.0, 6.0, 6.0, 3.0, 1.0, 2.0, 2.0, 7.0, 8.0, 2.0, 7.0, 10.0, 4.0, 8.0, 6.0, 7.0, 5.0, 7.0, 4.0, 9.0, 5.0, 3.0, 9.0, 6.0, 3.0, 6.0, 6.0, 9.0, 7.0, 4.0, 1.0, 5.0, 8.0, 10.0, 3.0, 4.0, 6.0, 6.0, 5.0, 5.0, 1.0, 4.0, 1.0, 7.0, 5.0, 9.0, 4.0, 8.0, 6.0, 6.0, 6.0, 8.0, 9.0, 7.0, 4.0, 4.0, 5.0, 6.0, 8.0, 8.0, 1.0, 7.0, 8.0, 7.0, 1.0, 4.0, 7.0, 9.0, 9.0, 8.0, 5.0, 4.0, 1.0, 1.0, 2.0, 6.0, 10.0, 10.0, 9.0, 4.0, 9.0, 10.0, 3.0, 4.0, 5.0, 9.0, 10.0, 3.0, 4.0, 5.0, 8.0, 1.0, 8.0, 1.0, 7.0, 1.0, 3.0, 1.0, 10.0, 4.0, 7.0, 6.0, 9.0, 7.0, 1.0, 5.0, 10.0, 6.0, 9.0, 10.0, 5.0, 4.0, 2.0, 2.0, 10.0, 4.0, 7.0, 7.0, 8.0, 8.0, 2.0, 5.0, 3.0, 1.0, 6.0, 3.0, 1.0, 2.0, 4.0, 6.0, 6.0]
global b_y = 10
global p = [0.712, 0.19, 0.65, 0.178, 0.885, 0.691, 0.715, 0.457, 0.909, 0.749, 0.251, 0.15, 0.032, 0.333, 0.604, 0.013, 0.769, 0.418, 0.546, 0.821, 0.627, 0.465, 0.308, 0.649, 0.657, 0.887, 0.543, 0.545, 0.582, 0.573, 0.317, 0.716, 0.198, 0.683, 0.35, 0.094, 0.99, 0.263, 0.613, 0.727, 0.027, 0.659, 0.651, 0.207, 0.18, 0.506, 0.053, 0.161, 0.877, 0.127, 0.699, 0.505, 0.461, 0.277, 0.221, 0.353, 0.27, 0.524, 0.288, 0.97, 0.968, 0.442, 0.268, 0.016, 0.061, 0.335, 0.439, 0.882, 0.717, 0.298, 0.046, 0.618, 0.22, 0.836, 0.043, 0.974, 0.134, 0.321, 0.996, 0.147, 0.127, 0.009, 0.342, 0.893, 0.727, 0.43, 0.364, 0.947, 0.487, 0.506, 0.025, 0.36, 0.548, 0.238, 0.715, 0.939, 0.95, 0.699, 0.479, 0.748, 0.071, 0.382, 0.944, 0.253, 0.547, 0.392, 0.643, 0.972, 0.29, 0.99, 0.162, 0.85, 0.758, 0.981, 0.354, 0.048, 0.878, 0.208, 0.972, 0.304, 0.506, 0.046, 0.862, 0.653, 0.018, 0.479, 0.786, 0.118, 0.535, 0.264, 0.431, 0.619, 0.966, 0.148, 0.393, 0.12, 0.101, 0.513, 0.354, 0.095, 0.465, 0.11, 0.08, 0.896, 0.297, 0.327, 0.985, 0.048, 0.395, 0.218, 0.627, 0.469, 0.725, 0.898, 0.423, 0.552, 0.56, 0.519, 0.497, 0.313, 0.405, 0.387, 0.951, 0.738, 0.24, 0.832, 0.632, 0.189, 0.196, 0.457, 0.461, 0.477, 0.138, 0.325, 0.719, 0.127, 0.17, 0.113, 0.788, 0.597, 0.521, 0.806, 0.25, 0.889, 0.935, 0.526, 0.811, 0.317, 0.058, 0.232, 0.178, 0.374, 0.26, 0.942, 0.55, 0.973, 0.065, 0.09, 0.953, 0.558, 0.838, 0.036, 0.022, 0.221, 0.354, 0.467, 0.959, 0.151, 0.471, 0.858, 0.086, 0.806, 0.311, 0.039, 0.131, 0.072, 0.991, 0.808, 0.564, 0.028, 0.015, 0.609, 0.302, 0.704, 0.625, 0.079, 0.417, 0.306, 0.562, 0.933, 0.011, 0.695, 0.977, 0.541, 0.81, 0.018, 0.293, 0.694, 0.37, 0.11, 0.864, 0.137, 0.677, 0.174, 0.502, 0.636, 0.821, 0.758, 0.926, 0.493, 0.71, 0.97, 0.508, 0.326, 0.996, 0.566, 0.326, 0.147, 0.038, 0.763, 0.609, 0.23]
global q = [0.857, 0.719, 0.961, 0.876, 0.973, 0.983, 0.977, 0.922, 0.928, 0.891, 0.33, 0.437, 0.076, 0.681, 0.836, 0.449, 0.842, 0.433, 0.873, 0.855, 0.897, 0.826, 0.603, 0.689, 0.723, 0.908, 0.912, 0.607, 0.758, 0.725, 0.965, 0.921, 0.258, 0.79, 0.695, 0.543, 0.997, 0.81, 0.962, 0.902, 0.179, 0.778, 0.686, 0.83, 0.827, 0.918, 0.655, 0.783, 0.99, 0.993, 0.752, 0.552, 0.612, 0.946, 0.488, 0.37, 0.547, 0.632, 0.9, 0.997, 0.979, 0.545, 0.941, 0.625, 0.614, 0.697, 0.912, 0.908, 0.992, 0.556, 0.567, 0.715, 0.903, 0.887, 0.894, 0.998, 0.158, 0.957, 0.999, 0.391, 0.968, 0.723, 0.554, 0.898, 0.95, 0.786, 0.826, 0.978, 0.568, 0.668, 0.398, 0.726, 0.744, 0.607, 0.935, 0.992, 0.989, 0.937, 0.601, 0.823, 0.744, 0.79, 0.957, 0.747, 0.873, 0.591, 0.755, 0.977, 0.712, 0.992, 0.21, 0.879, 0.919, 0.997, 0.988, 0.944, 0.901, 0.88, 0.979, 0.667, 0.88, 0.058, 0.943, 0.827, 0.698, 0.981, 0.906, 0.691, 0.657, 0.321, 0.633, 0.793, 0.972, 0.799, 0.473, 0.658, 0.485, 0.672, 0.473, 0.815, 0.562, 0.217, 0.747, 0.938, 0.693, 0.844, 0.99, 0.626, 0.675, 0.634, 0.67, 0.947, 0.858, 0.967, 0.466, 0.913, 0.759, 0.542, 0.606, 0.541, 0.769, 0.811, 0.993, 0.802, 0.684, 0.902, 0.794, 0.387, 0.589, 0.657, 0.52, 0.517, 0.824, 0.688, 0.761, 0.664, 0.74, 0.382, 0.795, 0.771, 0.62, 0.856, 0.271, 0.927, 0.995, 0.814, 0.909, 0.349, 0.965, 0.893, 0.956, 0.425, 0.273, 0.998, 0.846, 0.98, 0.171, 0.48, 0.975, 0.899, 0.995, 0.57, 0.984, 0.489, 0.923, 0.939, 0.965, 0.182, 0.891, 0.922, 0.665, 0.827, 0.759, 0.237, 0.569, 0.084, 0.997, 0.82, 0.91, 0.167, 0.3, 0.733, 0.441, 0.902, 0.651, 0.597, 0.458, 0.998, 0.638, 0.977, 0.98, 0.915, 0.991, 0.645, 0.943, 0.373, 0.764, 0.906, 0.754, 0.367, 0.897, 0.518, 0.953, 0.851, 0.572, 0.79, 0.871, 0.929, 0.942, 0.788, 0.777, 0.997, 0.886, 0.934, 0.996, 0.84, 0.404, 0.398, 0.253, 0.859, 0.829, 0.722]
global origin = 1
global destination = 50