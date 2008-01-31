#     This file is part of FAiR, a program to conduct Factor Analysis in R
#     Copyright 2008 Benjamin King Goodrich
# 
#     FAiR is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     FAiR is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with FAiR.  If not, see <http://www.gnu.org/licenses/>.

.onAttach <- function( ... ) {
FAiRLib <- dirname(system.file(package = "FAiR"))
version <- packageDescription("FAiR", lib = FAiRLib)$Version
BuildDate <- packageDescription("FAiR", lib = FAiRLib)$Date
cat(paste("\n##  FAiR Version", version, "Build Date:", BuildDate, "\n"))
cat("##  See http://wiki.r-project.org/rwiki/doku.php?id=packages:cran:fair for more information.\n")
cat("FAiR  Copyright (C) 2008  Benjamin King Goodrich\n")
cat("This program comes with ABSOLUTELY NO WARRANTY.\n")
cat("This is free software, and you are welcome to redistribute it\n")
cat("under certain conditions, namely those specified at\n")
cat("http://www.gnu.org/licenses/gpl.txt\n")

if (length(grep("darwin", R.version$platform))) {
      cat("\n\nWARNING: It appears you are using a Mac.\n",
      "FAiR will CRASH under normal usage unless the X server is running.\n",
      "Verify that the X server before calling any function in FAiR that uses a GUI menu.\n",
      "To do so, click the big X icon in the middle of the R GUI or\n",
      "start the R (GUI) session from xterm after starting the X server manually\n",
      "\n\n")
}

else if (.Platform$OS.type == "windows") {
	flush.console()
	cat("\n\nIt appears you are using Windows.\n",
	    "It is recommended that you disable buffering by pressing Ctrl-W or\n",
	    "by deselecting Misc -> Buffered output in the menu at the top.\n",
	    "Doing so will permit more consistent printing of the progress of the genetic algorithm.\n")
	flush.console()
}
}

