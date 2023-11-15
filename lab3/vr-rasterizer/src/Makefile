all: fullscreen.vert.spv debug_pattern.frag.spv
fullscreen.vert.spv: fullscreen.vert
	/usr/bin/glslc fullscreen.vert -o fullscreen.vert.spv

debug_pattern.frag.spv: debug_pattern.frag
	/usr/bin/glslc debug_pattern.frag -o debug_pattern.frag.spv

clean:
	rm -rf *.spv