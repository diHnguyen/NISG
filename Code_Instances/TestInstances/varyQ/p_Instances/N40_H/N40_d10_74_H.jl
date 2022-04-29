global arcs = [1 3; 1 4; 1 20; 2 6; 2 7; 2 11; 2 38; 3 4; 3 7; 3 16; 3 23; 4 23; 5 6; 5 9; 5 19; 6 12; 7 10; 8 2; 8 6; 8 29; 8 32; 9 16; 9 19; 9 21; 9 25; 9 29; 10 31; 10 36; 10 40; 11 13; 11 21; 11 29; 12 21; 12 28; 13 3; 13 7; 13 18; 13 30; 13 38; 14 7; 14 10; 14 22; 14 24; 15 11; 15 16; 15 17; 15 22; 15 35; 16 6; 16 14; 16 33; 17 10; 17 14; 17 23; 18 3; 18 8; 18 11; 18 15; 18 16; 18 25; 18 36; 18 37; 18 39; 19 30; 19 34; 20 5; 20 13; 20 22; 20 27; 20 34; 21 9; 22 14; 22 27; 22 39; 23 5; 23 7; 23 8; 23 9; 24 18; 24 30; 25 17; 26 7; 26 27; 26 29; 27 7; 28 6; 28 13; 28 14; 28 31; 28 39; 29 34; 30 2; 30 21; 31 3; 31 4; 31 18; 31 36; 32 13; 32 14; 32 21; 32 40; 33 2; 33 14; 33 17; 33 26; 33 36; 33 38; 33 40; 34 19; 34 20; 34 24; 34 29; 35 18; 35 22; 35 29; 36 2; 36 15; 36 16; 36 28; 36 32; 37 10; 37 23; 37 35; 38 8; 38 37; 39 3; 39 9; 39 11; 39 14; 39 16; 39 26; 39 40]
global d_x = [8.0, 3.0, 5.0, 1.0, 10.0, 3.0, 8.0, 8.0, 3.0, 3.0, 6.0, 1.0, 10.0, 1.0, 6.0, 9.0, 7.0, 2.0, 6.0, 4.0, 10.0, 6.0, 6.0, 4.0, 9.0, 7.0, 8.0, 2.0, 1.0, 9.0, 7.0, 1.0, 4.0, 6.0, 10.0, 2.0, 9.0, 4.0, 4.0, 7.0, 10.0, 9.0, 6.0, 2.0, 6.0, 4.0, 8.0, 3.0, 8.0, 3.0, 4.0, 1.0, 1.0, 1.0, 1.0, 6.0, 6.0, 8.0, 10.0, 8.0, 6.0, 3.0, 3.0, 4.0, 3.0, 2.0, 7.0, 7.0, 4.0, 8.0, 3.0, 1.0, 10.0, 2.0, 9.0, 8.0, 2.0, 7.0, 9.0, 10.0, 1.0, 7.0, 10.0, 2.0, 9.0, 10.0, 2.0, 2.0, 5.0, 7.0, 2.0, 1.0, 9.0, 7.0, 2.0, 4.0, 7.0, 10.0, 2.0, 2.0, 6.0, 4.0, 6.0, 8.0, 6.0, 2.0, 9.0, 3.0, 6.0, 9.0, 2.0, 10.0, 6.0, 8.0, 2.0, 3.0, 7.0, 7.0, 9.0, 2.0, 4.0, 5.0, 4.0, 7.0, 4.0, 10.0, 7.0, 5.0, 8.0, 3.0, 2.0, 4.0]
global b_x = 5
global d_y = [3.0, 7.0, 4.0, 1.0, 1.0, 10.0, 8.0, 1.0, 5.0, 9.0, 3.0, 3.0, 7.0, 9.0, 7.0, 4.0, 3.0, 10.0, 2.0, 2.0, 1.0, 7.0, 1.0, 9.0, 8.0, 4.0, 6.0, 5.0, 6.0, 8.0, 3.0, 4.0, 3.0, 2.0, 9.0, 9.0, 4.0, 10.0, 7.0, 6.0, 5.0, 3.0, 8.0, 6.0, 3.0, 1.0, 9.0, 4.0, 9.0, 6.0, 3.0, 1.0, 5.0, 3.0, 8.0, 3.0, 4.0, 4.0, 10.0, 10.0, 10.0, 8.0, 6.0, 8.0, 6.0, 3.0, 1.0, 6.0, 5.0, 3.0, 5.0, 4.0, 7.0, 1.0, 4.0, 1.0, 2.0, 7.0, 5.0, 2.0, 7.0, 5.0, 1.0, 8.0, 2.0, 4.0, 6.0, 8.0, 1.0, 1.0, 3.0, 8.0, 5.0, 9.0, 1.0, 2.0, 6.0, 2.0, 3.0, 9.0, 2.0, 6.0, 6.0, 1.0, 7.0, 7.0, 10.0, 6.0, 9.0, 3.0, 8.0, 3.0, 1.0, 4.0, 5.0, 4.0, 9.0, 3.0, 7.0, 4.0, 10.0, 9.0, 5.0, 7.0, 2.0, 3.0, 7.0, 3.0, 8.0, 2.0, 8.0, 6.0]
global b_y = 10
global p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
global q = [0.861, 0.948, 0.892, 0.926, 0.828, 0.911, 0.947, 0.918, 0.915, 0.975, 0.987, 0.812, 0.887, 0.892, 0.992, 0.985, 0.94, 0.963, 0.972, 0.996, 0.86, 0.877, 0.802, 0.809, 0.916, 0.842, 0.884, 0.809, 0.897, 0.885, 0.898, 0.937, 0.943, 0.928, 0.86, 0.86, 0.841, 0.87, 0.881, 0.878, 0.989, 0.943, 0.833, 0.889, 0.835, 0.982, 0.86, 0.982, 0.947, 0.847, 0.997, 0.88, 0.933, 0.969, 0.864, 0.966, 0.964, 0.817, 0.937, 0.92, 0.973, 0.845, 0.912, 0.997, 0.877, 0.965, 0.986, 0.923, 0.853, 0.845, 0.807, 0.846, 0.911, 0.892, 0.808, 0.854, 0.876, 0.806, 0.875, 0.884, 0.92, 0.8, 0.848, 0.808, 0.877, 0.99, 0.87, 0.855, 0.999, 0.897, 0.853, 0.809, 0.869, 0.929, 0.968, 0.836, 0.903, 0.86, 0.885, 0.838, 0.825, 0.914, 0.925, 0.882, 0.942, 0.97, 0.894, 0.821, 0.865, 0.952, 0.955, 0.82, 0.94, 0.846, 0.847, 0.935, 0.807, 0.898, 0.973, 0.873, 0.984, 0.885, 0.914, 0.999, 0.835, 0.898, 0.831, 0.875, 0.837, 0.876, 0.864, 0.882]
global origin = 1
global destination = 40