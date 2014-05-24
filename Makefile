PREFIX=/usr/local

all:
	@bash checkinstall.sh
	@mkdir -p build
	@cd build && cmake -D "CMAKE_INSTALL_PREFIX=$(PREFIX)" ../
	@echo ""
	@echo ""
	@echo "To compile and install to $(PREFIX), type: "
	@echo ""
	@echo "cd build && make install"
	@echo ""
	

clean:
	@rm -rf build/
