OBJS = camd_i_aat.o camd_i_1.o camd_i_2.o camd_i_dump.o camd_i_postorder.o camd_i_defaults.o camd_i_order.o camd_i_control.o camd_i_info.o camd_i_valid.o camd_i_preprocess.o camd_l_aat.o camd_l_1.o camd_l_2.o camd_l_dump.o camd_l_postorder.o camd_l_defaults.o camd_l_order.o camd_l_control.o camd_l_info.o camd_l_valid.o camd_l_preprocess.o camd_global.o

camd_i_aat.o: camd_aat.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDINT -c camd_aat.c -o $@
camd_i_1.o: camd_1.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDINT -c camd_1.c -o $@
camd_i_2.o: camd_2.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDINT -c camd_2.c -o $@
camd_i_dump.o: camd_dump.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDINT -c camd_dump.c -o $@
camd_i_postorder.o: camd_postorder.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDINT -c camd_postorder.c -o $@
camd_i_defaults.o: camd_defaults.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDINT -c camd_defaults.c -o $@
camd_i_order.o: camd_order.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDINT -c camd_order.c -o $@
camd_i_control.o: camd_control.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDINT -c camd_control.c -o $@
camd_i_info.o: camd_info.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDINT -c camd_info.c -o $@
camd_i_valid.o: camd_valid.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDINT -c camd_valid.c -o $@
camd_i_preprocess.o: camd_preprocess.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDINT -c camd_preprocess.c -o $@


camd_l_aat.o: camd_aat.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDLONG -c camd_aat.c -o $@
camd_l_1.o: camd_1.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDLONG -c camd_1.c -o $@
camd_l_2.o: camd_2.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDLONG -c camd_2.c -o $@
camd_l_dump.o: camd_dump.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDLONG -c camd_dump.c -o $@
camd_l_postorder.o: camd_postorder.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDLONG -c camd_postorder.c -o $@
camd_l_defaults.o: camd_defaults.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDLONG -c camd_defaults.c -o $@
camd_l_order.o: camd_order.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDLONG -c camd_order.c -o $@
camd_l_control.o: camd_control.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDLONG -c camd_control.c -o $@
camd_l_info.o: camd_info.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDLONG -c camd_info.c -o $@
camd_l_valid.o: camd_valid.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDLONG -c camd_valid.c -o $@
camd_l_preprocess.o: camd_preprocess.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDLONG -c camd_preprocess.c -o $@
camd_global.o: camd_global.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -c camd_global.c -o $@

