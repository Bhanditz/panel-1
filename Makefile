pkg = panel

all: Makefile DESCRIPTION
	@echo "building package \`$(pkg)'"


mostlyclean: clean
clean:
	@if test -d src; then (cd src && $(MAKE) $@); fi
	@if test -d doc; then (cd doc && $(MAKE) $@); fi
