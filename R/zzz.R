.onLoad <- function(lib, pkg) {
if(.Platform$OS.type == "windows" && interactive()  && .Platform$GUI == "Rgui") {
	addVigs2WinMenu("GEOsubmission")
	}

}