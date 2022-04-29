global arcs = [1 4; 1 13; 1 26; 2 9; 2 10; 2 13; 2 16; 2 19; 2 25; 2 28; 2 30; 3 5; 3 13; 3 17; 4 14; 4 20; 5 8; 6 7; 6 14; 7 6; 7 8; 7 15; 8 3; 8 4; 8 14; 8 21; 9 11; 9 22; 9 29; 10 4; 10 23; 11 2; 11 3; 11 6; 11 12; 11 18; 12 4; 12 6; 12 14; 12 16; 12 19; 13 4; 13 9; 13 11; 13 27; 13 29; 14 9; 14 12; 15 2; 15 6; 16 28; 17 9; 17 10; 17 13; 17 22; 17 24; 17 29; 18 3; 18 6; 18 9; 18 10; 19 6; 19 18; 19 21; 19 27; 20 9; 20 27; 21 22; 22 6; 23 5; 23 15; 24 26; 25 7; 26 3; 26 27; 27 11; 27 12; 27 15; 27 21; 28 29; 29 24]
global d_x = [5.0, 10.0, 10.0, 1.0, 2.0, 7.0, 6.0, 6.0, 4.0, 2.0, 8.0, 8.0, 9.0, 2.0, 1.0, 7.0, 8.0, 1.0, 6.0, 9.0, 1.0, 8.0, 7.0, 7.0, 4.0, 5.0, 3.0, 6.0, 6.0, 5.0, 2.0, 8.0, 3.0, 4.0, 5.0, 2.0, 6.0, 5.0, 2.0, 1.0, 4.0, 10.0, 10.0, 3.0, 9.0, 4.0, 4.0, 3.0, 6.0, 2.0, 9.0, 2.0, 3.0, 1.0, 8.0, 3.0, 5.0, 1.0, 1.0, 9.0, 10.0, 9.0, 3.0, 7.0, 1.0, 5.0, 5.0, 7.0, 2.0, 8.0, 1.0, 5.0, 5.0, 9.0, 2.0, 1.0, 5.0, 7.0, 4.0, 9.0, 3.0]
global b_x = 5
global d_y = [9.0, 2.0, 4.0, 3.0, 3.0, 4.0, 8.0, 8.0, 3.0, 1.0, 10.0, 5.0, 10.0, 3.0, 7.0, 5.0, 10.0, 1.0, 6.0, 4.0, 3.0, 4.0, 2.0, 8.0, 9.0, 10.0, 5.0, 10.0, 7.0, 4.0, 3.0, 9.0, 2.0, 4.0, 6.0, 6.0, 1.0, 10.0, 1.0, 2.0, 8.0, 2.0, 4.0, 5.0, 8.0, 7.0, 2.0, 8.0, 7.0, 8.0, 9.0, 8.0, 10.0, 1.0, 6.0, 1.0, 7.0, 10.0, 6.0, 10.0, 6.0, 2.0, 1.0, 9.0, 2.0, 8.0, 1.0, 3.0, 9.0, 10.0, 4.0, 8.0, 4.0, 4.0, 10.0, 9.0, 1.0, 4.0, 2.0, 10.0, 6.0]
global b_y = 10
global p = [0.662, 0.65, 0.312, 0.138, 0.167, 0.972, 0.65, 0.906, 0.745, 0.304, 0.668, 0.811, 0.866, 0.642, 0.459, 0.908, 0.425, 0.216, 0.072, 0.226, 0.458, 0.915, 0.605, 0.538, 0.049, 0.51, 0.316, 0.927, 0.021, 0.624, 0.048, 0.718, 0.747, 0.443, 0.793, 0.522, 0.309, 0.865, 0.001, 0.81, 0.191, 0.209, 0.974, 0.11, 0.938, 0.424, 0.843, 0.58, 0.71, 0.616, 0.279, 0.662, 0.409, 0.862, 0.815, 0.363, 0.749, 0.999, 0.534, 0.421, 0.846, 0.281, 0.83, 0.144, 0.093, 0.954, 0.031, 0.916, 0.435, 0.855, 0.309, 0.53, 0.324, 0.14, 0.866, 0.634, 0.399, 0.603, 0.787, 0.785, 0.237]
global q = [0.767, 0.939, 0.344, 0.793, 0.31, 0.99, 0.868, 0.975, 0.858, 0.501, 0.761, 0.95, 0.874, 0.687, 0.569, 0.974, 0.469, 0.624, 0.127, 0.661, 0.735, 0.965, 0.796, 0.576, 0.182, 0.939, 0.509, 0.951, 0.174, 0.885, 0.559, 0.988, 0.841, 0.517, 0.877, 0.831, 0.316, 0.951, 0.624, 0.905, 0.858, 0.442, 0.978, 0.221, 0.969, 0.69, 0.966, 0.788, 0.952, 0.714, 0.418, 0.738, 0.472, 0.865, 0.925, 0.979, 0.859, 0.999, 0.837, 0.552, 0.877, 0.907, 0.942, 0.833, 0.171, 0.985, 0.386, 0.922, 0.809, 0.906, 0.724, 0.559, 0.639, 0.231, 0.96, 0.725, 0.613, 0.73, 0.907, 0.944, 0.237]
global origin = 1
global destination = 30
