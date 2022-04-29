global arcs = [1 7; 1 12; 1 30; 1 40; 2 9; 2 34; 3 2; 3 23; 3 31; 3 32; 3 35; 4 16; 4 29; 5 9; 5 21; 5 28; 5 36; 6 10; 6 12; 6 23; 6 31; 6 33; 7 11; 7 25; 7 27; 7 31; 7 33; 7 35; 8 10; 8 18; 8 28; 8 39; 9 4; 9 12; 9 18; 9 34; 10 8; 10 31; 11 8; 11 31; 11 36; 11 38; 12 6; 12 7; 12 8; 12 24; 12 27; 12 33; 13 14; 13 15; 13 26; 13 40; 14 8; 14 11; 14 25; 14 27; 15 3; 15 9; 15 12; 16 18; 16 19; 17 4; 17 8; 17 24; 17 32; 17 33; 17 35; 18 3; 18 8; 18 37; 19 10; 19 22; 19 39; 20 7; 20 8; 20 9; 20 23; 21 2; 21 6; 21 13; 21 14; 21 17; 21 29; 21 32; 22 10; 23 25; 23 34; 24 11; 24 12; 24 13; 24 23; 24 34; 24 35; 25 7; 25 21; 25 28; 26 18; 26 31; 26 35; 26 36; 26 38; 27 11; 27 12; 27 15; 27 17; 27 22; 27 28; 28 7; 28 31; 28 33; 28 34; 28 35; 29 9; 30 3; 30 15; 30 26; 30 28; 30 40; 31 3; 31 8; 31 11; 31 18; 31 19; 31 35; 32 37; 32 38; 33 3; 33 4; 33 5; 33 6; 33 16; 33 17; 33 18; 33 30; 34 33; 35 12; 35 20; 35 22; 35 25; 35 33; 36 13; 36 29; 37 4; 37 13; 38 4; 38 5; 38 7; 38 17; 39 15; 39 25]
global d_x = [3.0, 6.0, 5.0, 5.0, 1.0, 6.0, 1.0, 6.0, 2.0, 3.0, 6.0, 4.0, 9.0, 1.0, 4.0, 5.0, 6.0, 8.0, 8.0, 3.0, 1.0, 9.0, 3.0, 5.0, 7.0, 6.0, 3.0, 8.0, 10.0, 10.0, 8.0, 3.0, 6.0, 4.0, 8.0, 9.0, 8.0, 7.0, 2.0, 3.0, 7.0, 8.0, 5.0, 6.0, 4.0, 4.0, 2.0, 2.0, 4.0, 2.0, 7.0, 4.0, 3.0, 9.0, 3.0, 1.0, 9.0, 10.0, 5.0, 4.0, 10.0, 6.0, 8.0, 3.0, 5.0, 4.0, 7.0, 5.0, 8.0, 9.0, 9.0, 7.0, 9.0, 4.0, 10.0, 2.0, 10.0, 6.0, 3.0, 7.0, 2.0, 1.0, 7.0, 5.0, 1.0, 6.0, 1.0, 6.0, 1.0, 7.0, 6.0, 6.0, 3.0, 5.0, 6.0, 2.0, 2.0, 4.0, 3.0, 3.0, 6.0, 6.0, 10.0, 1.0, 2.0, 6.0, 1.0, 3.0, 10.0, 3.0, 10.0, 3.0, 10.0, 4.0, 1.0, 7.0, 8.0, 10.0, 10.0, 5.0, 10.0, 7.0, 9.0, 2.0, 5.0, 8.0, 8.0, 5.0, 10.0, 2.0, 4.0, 2.0, 5.0, 6.0, 10.0, 9.0, 2.0, 2.0, 6.0, 2.0, 5.0, 1.0, 5.0, 3.0, 4.0, 4.0, 1.0, 6.0, 9.0, 10.0]
global b_x = 5
global d_y = [8.0, 4.0, 8.0, 1.0, 3.0, 10.0, 5.0, 9.0, 5.0, 9.0, 5.0, 5.0, 1.0, 7.0, 3.0, 1.0, 10.0, 6.0, 6.0, 9.0, 7.0, 1.0, 4.0, 3.0, 6.0, 2.0, 6.0, 9.0, 9.0, 9.0, 1.0, 6.0, 10.0, 10.0, 1.0, 10.0, 7.0, 4.0, 8.0, 3.0, 10.0, 6.0, 1.0, 10.0, 8.0, 8.0, 6.0, 8.0, 5.0, 2.0, 5.0, 10.0, 7.0, 3.0, 8.0, 7.0, 6.0, 10.0, 5.0, 5.0, 1.0, 6.0, 8.0, 1.0, 4.0, 5.0, 1.0, 9.0, 7.0, 7.0, 3.0, 2.0, 7.0, 4.0, 10.0, 8.0, 2.0, 1.0, 8.0, 5.0, 7.0, 5.0, 3.0, 2.0, 8.0, 5.0, 6.0, 4.0, 5.0, 7.0, 2.0, 7.0, 4.0, 9.0, 9.0, 4.0, 7.0, 2.0, 4.0, 3.0, 1.0, 6.0, 8.0, 1.0, 3.0, 3.0, 7.0, 5.0, 3.0, 5.0, 9.0, 8.0, 10.0, 8.0, 8.0, 10.0, 7.0, 5.0, 9.0, 9.0, 3.0, 10.0, 2.0, 5.0, 7.0, 5.0, 8.0, 1.0, 4.0, 7.0, 2.0, 7.0, 4.0, 3.0, 8.0, 2.0, 9.0, 4.0, 1.0, 1.0, 5.0, 8.0, 8.0, 4.0, 8.0, 8.0, 2.0, 9.0, 3.0, 9.0]
global b_y = 10
global p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
global q = [0.994, 0.978, 0.811, 0.945, 0.96, 0.838, 0.872, 0.896, 0.966, 0.879, 0.954, 0.838, 0.856, 0.86, 0.879, 0.903, 0.805, 0.837, 0.901, 0.829, 0.957, 0.865, 0.87, 0.912, 0.812, 0.869, 0.814, 0.856, 0.858, 0.993, 0.843, 0.892, 0.998, 0.914, 0.832, 0.857, 0.869, 0.86, 0.948, 0.851, 0.862, 0.995, 0.878, 0.863, 0.883, 0.911, 0.973, 0.841, 0.841, 0.856, 0.986, 0.954, 0.99, 0.853, 0.843, 0.806, 0.855, 0.884, 0.984, 0.852, 0.959, 0.954, 0.92, 0.87, 0.804, 0.86, 0.932, 0.991, 0.867, 0.81, 0.985, 0.908, 0.997, 0.94, 0.982, 0.875, 0.92, 0.833, 0.854, 0.924, 0.971, 0.987, 0.848, 0.919, 0.841, 0.959, 0.808, 0.914, 0.852, 0.842, 0.851, 0.864, 0.919, 0.855, 0.897, 0.846, 0.908, 0.942, 0.965, 0.921, 0.835, 0.844, 0.953, 0.938, 0.975, 0.989, 0.88, 0.846, 0.96, 0.891, 0.864, 0.918, 0.819, 0.954, 0.837, 0.891, 0.846, 0.912, 0.929, 0.979, 0.972, 0.833, 0.832, 0.935, 0.924, 0.836, 0.933, 0.891, 0.839, 0.867, 0.851, 0.949, 0.857, 0.976, 0.929, 0.906, 0.806, 0.996, 0.845, 0.822, 0.834, 0.833, 0.943, 0.876, 0.945, 0.974, 0.807, 0.8, 0.857, 0.994]
global origin = 1
global destination = 40