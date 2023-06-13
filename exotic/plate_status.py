
class PlateStatus:
    def __init__(self, logfunc):
        self.statusByFilename = dict()
        self.filenameList = []
        self.filename = "N/A"
        self.errorcodes = set()
        self.logfunc = logfunc
        self.errorcodes.add("outofframe_target")
        self.errorcodes.add("lowflux_target")
        self.errorcodes.add("skybg_target")
        self.errorcodes.add("fits_error")
        self.errorcodes.add("alignment_error")
    # Initialze set of files, as well as ordered index
    def initializeFilenames(self, filenames: list):
        self.filenameList = filenames.copy()
        self.filenameList.sort()
        for fneme in filenames:
            if fneme not in self.statusByFilename:
                self.statusByFilename[fneme] = {}
        return self
    # Initialize comparison star count
    def initializeComparisonStarCount(self, compCount: int):
        for i in range(compCount):
            self.errorcodes.add(f"outofframe_comp{i+1}")
            self.errorcodes.add(f"lowflux_comp{i+1}")
            self.errorcodes.add(f"skybg_comp{i+1}")            
    # Sets current filename (for any reported errors) - sets starIndex=0 (target)
    def setCurrentFilename(self, filename: str):
        filename = str(filename)
        if filename not in self.statusByFilename:
            self.statusByFilename[filename] = {}
        self.filename = filename
        return self
    # Log an error
    def _logError(self, errorcode: str, message: str) -> None:
        if self.filename not in self.statusByFilename:
            self.statusByFilename[self.filename] = {}
        rec = self.statusByFilename[self.filename]
        if errorcode in rec:
            return
        # Mark error on this file
        rec[errorcode] = True
        self.errorcodes.add(errorcode)
        # And log new warning
        self.logfunc(message, warn=True)
    # Report out of frame warning for start ;index' (0=target, 1+=comp #N)
    def outOfFrameWarning(self, starIndex):
        if starIndex == 0:  # Target star
            self._logError("outofframe_target",
                f"Target star beyond edge of file {self.filename}")
        else:
            self._logError(f"outofframe_comp{starIndex}",
                f"Comparison star #{starIndex} star beyond edge of file {self.filename}")
    # Report low flux amplitude warning for start ;index' (0=target, 1+=comp #N)
    def lowFluxAmplitudeWarning(self, starIndex: int, xc: float, yc: float):
        if starIndex == 0:  # Target star
            self._logError("lowflux_target",
                f"Measured flux for Target star is low in file  {self.filename} - are you sure there is a star at [{xc:.1f}, {yc:.1f}]?")
        else:
            self._logError(f"lowflux_comp{starIndex}",
                f"Measured flux for Comparison star #{starIndex} is low in file {self.filename} - are you sure there is a star at [{xc:.1f}, {yc:.1f}]?")
    # Report sky background warning for start ;index' (0=target, 1+=comp #N)
    def skyBackgroundWarning(self, starIndex: int, xc: float, yc: float):
        if starIndex == 0:  # Target star
            self._logError("skybg_target",
                f"Sky background error for Target star for file {self.filename} - are you sure there is a star at [{xc:.1f}, {yc:.1f}]?")
        else:
            self._logError(f"skybg_comp{starIndex}",
                f"Sky background error for Comparison star #{starIndex} for file {self.filename} - are you sure there is a star at [{xc:.1f}, {yc:.1f}]?")
    # Reort file format error
    def fitsFormatError(self, e: OSError):
        self._logError("fits_error",
            f"Corrupted file {self.filename} ({e}) --removed from reduction")
    # Reort alignment error
    def setObsTime(self, time):
        if self.filename not in self.statusByFilename:
            self.statusByFilename[self.filename] = {}
        self.statusByFilename[self.filename]['time'] = time
        return self
    # Reort frame time
    def alignmentError(self):
        self._logError("alignment_error",
            f"File {self.filename} failed to align with first file")
    # Write plate status to CSV file
    def writePlateStatus(self, file: str):
        with open(file, 'w') as f:
            cols = list(self.errorcodes)
            cols.sort()
            f.write(f"# filename,time,{','.join(cols)}\n")
            for file in self.filenameList:
                rec = self.statusByFilename[file]
                line = f"\"{file}\",{rec['time'] if 'time' in rec else ''}"
                for col in cols:
                    if col in rec:
                        line = f"{line},{rec[col]}"
                    else:
                        line = f"{line},False"
                f.write(f"{line}\n")
