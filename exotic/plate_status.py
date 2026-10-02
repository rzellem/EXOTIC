
class PlateStatus:
    WARNING_PROGRESS_INTERVAL = 100
    WARNING_CONDITION_TEXT = {
        "outofframe": "outside the image",
        "lowflux": "low flux",
        "overexposed": "overexposed",
        "skybg": "a sky-background failure",
    }

    def __init__(self, logfunc):
        self.statusByFilename = dict()
        self.filenameList = []
        self.filename = "N/A"
        self.errorcodes = set()
        self.logfunc = logfunc
        self.comparisonStarLabels = {}
        self.aggregatedWarnings = {}
        self.lastAggregatedWarningSummary = {}
        self.aggregationNoticeLogged = False
        self.errorcodes.add("outofframe_target")
        self.errorcodes.add("lowflux_target")
        self.errorcodes.add("overexposed_target")
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
            self.errorcodes.add(f"overexposed_comp{i+1}")
            self.errorcodes.add(f"skybg_comp{i+1}")

    def setComparisonStarLabels(self, labels=None):
        normalized_labels = {}
        for star_index, label in (labels or {}).items():
            try:
                normalized_index = int(star_index)
            except (TypeError, ValueError):
                continue
            if normalized_index <= 0 or not isinstance(label, str) or not label.strip():
                continue
            normalized_labels[normalized_index] = label.strip()
        self.comparisonStarLabels = normalized_labels
        return self

    def _starLabel(self, starIndex: int):
        if starIndex == 0:
            return "Target star"
        return self.comparisonStarLabels.get(starIndex, f"Comparison star #{starIndex}")

    def _displayFilename(self):
        return str(self.filename).replace('\\', '/').rsplit('/', 1)[-1]

    def _conditionCountText(self, condition: str, count: int):
        condition_text = self.WARNING_CONDITION_TEXT.get(condition, condition)
        return f"{condition_text} in {count} frame(s)"

    # Sets current filename (for any reported errors) - sets starIndex=0 (target)
    def setCurrentFilename(self, filename: str):
        filename = str(filename)
        if filename not in self.statusByFilename:
            self.statusByFilename[filename] = {}
        self.filename = filename
        return self
    # Log an error
    def _logError(self, errorcode: str, message: str, starIndex: int = None,
                  condition: str = None, starLabel: str = None) -> None:
        if self.filename not in self.statusByFilename:
            self.statusByFilename[self.filename] = {}
        rec = self.statusByFilename[self.filename]
        if errorcode in rec:
            return
        # Mark error on this file
        rec[errorcode] = True
        self.errorcodes.add(errorcode)
        if starIndex is None or condition is None:
            self.logfunc(message, warn=True)
            return

        label = (
            starLabel.strip()
            if isinstance(starLabel, str) and starLabel.strip()
            else self._starLabel(starIndex)
        )
        aggregate = self.aggregatedWarnings.setdefault(errorcode, {
            'star_index': starIndex,
            'label': label,
            'condition': condition,
            'count': 0,
        })
        aggregate['label'] = label
        aggregate['count'] += 1
        warning_count = aggregate['count']

        if warning_count == 1:
            self.logfunc(message, warn=True)
            if not self.aggregationNoticeLogged:
                self.logfunc(
                    "Plate-status detail: repeated frame-level star warnings are aggregated after "
                    "their first occurrence; running counts are reported every "
                    f"{self.WARNING_PROGRESS_INTERVAL} frames and exact per-frame flags are preserved "
                    "in the PlateStatus CSV."
                )
                self.aggregationNoticeLogged = True
        elif warning_count % self.WARNING_PROGRESS_INTERVAL == 0:
            self.logfunc(
                "Plate-status warning update: "
                f"{label}: {self._conditionCountText(condition, warning_count)} so far.",
                warn=True,
            )

    def logAggregatedWarningSummary(self):
        current_counts = {
            errorcode: aggregate['count']
            for errorcode, aggregate in self.aggregatedWarnings.items()
        }
        if current_counts == self.lastAggregatedWarningSummary:
            return

        repeated = [
            aggregate for aggregate in self.aggregatedWarnings.values()
            if aggregate['count'] > 1
        ]
        self.lastAggregatedWarningSummary = current_counts
        if not repeated:
            return

        grouped = {}
        for aggregate in repeated:
            group = grouped.setdefault(aggregate['star_index'], {
                'label': aggregate['label'],
                'conditions': [],
                'total': 0,
            })
            group['label'] = aggregate['label']
            group['conditions'].append(
                self._conditionCountText(aggregate['condition'], aggregate['count'])
            )
            group['total'] += aggregate['count']

        self.logfunc(
            "Plate-status warning summary: repeated frame-level diagnostics were aggregated; "
            "exact per-frame flags are preserved in the PlateStatus CSV."
        )
        for group in sorted(
            grouped.values(),
            key=lambda item: (-item['total'], item['label']),
        ):
            self.logfunc(f">-- {group['label']}: {'; '.join(group['conditions'])}.")

    # Report out of frame warning for start ;index' (0=target, 1+=comp #N)
    def outOfFrameWarning(self, starIndex):
        label = self._starLabel(starIndex)
        errorcode = "outofframe_target" if starIndex == 0 else f"outofframe_comp{starIndex}"
        self._logError(
            errorcode,
            f"{label} is beyond the edge of file {self._displayFilename()}",
            starIndex=starIndex,
            condition="outofframe",
            starLabel=label,
        )
    # Report low flux amplitude warning for start ;index' (0=target, 1+=comp #N)
    def lowFluxAmplitudeWarning(self, starIndex: int, xc: float, yc: float):
        label = self._starLabel(starIndex)
        errorcode = "lowflux_target" if starIndex == 0 else f"lowflux_comp{starIndex}"
        self._logError(
            errorcode,
            f"Measured flux for {label} is low in file {self._displayFilename()} - "
            f"are you sure there is a star at [{xc:.1f}, {yc:.1f}]?",
            starIndex=starIndex,
            condition="lowflux",
            starLabel=label,
        )
    # Report overexposure warning for star index (0=target, 1+=comp #N)
    def overexposedWarning(self, starIndex: int, xc: float, yc: float, threshold: float,
                           starLabel: str = None):
        label = (
            starLabel.strip()
            if isinstance(starLabel, str) and starLabel.strip()
            else self._starLabel(starIndex)
        )
        errorcode = "overexposed_target" if starIndex == 0 else f"overexposed_comp{starIndex}"
        self._logError(
            errorcode,
            f"{label} is overexposed in file {self._displayFilename()}; aperture pixels near "
            f"[{xc:.1f}, {yc:.1f}] exceeded {threshold:.1f}.",
            starIndex=starIndex,
            condition="overexposed",
            starLabel=label,
        )
    # Report sky background warning for start ;index' (0=target, 1+=comp #N)
    def skyBackgroundWarning(self, starIndex: int, xc: float, yc: float):
        label = self._starLabel(starIndex)
        errorcode = "skybg_target" if starIndex == 0 else f"skybg_comp{starIndex}"
        self._logError(
            errorcode,
            f"Sky background error for {label} in file {self._displayFilename()} - "
            f"are you sure there is a star at [{xc:.1f}, {yc:.1f}]?",
            starIndex=starIndex,
            condition="skybg",
            starLabel=label,
        )
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
        self.logAggregatedWarningSummary()
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
