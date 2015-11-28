include $(MkInclude)

PKG_CPPFLAGS = -I../Include -I../../SuiteSparse_config

LIB = ../../CAMD.a

lib: $(LIB)

include make_o.mk

$(LIB): $(OBJS)
	$(AR) -rcs $(LIB) $(OBJS)

mostlyclean: clean

clean:
	@-rm -rf .libs _libs $(LIB)
	@-rm -f *.o