// Copyright 2015 Ben H. Roos <beenhroos@gmail.com>
// Copyright 2016 Reed A. Cartwright <reed@cartwrig.ht>
// Copyright 2026 Robert R. Fitak <robert.fitak@ucf.edu>

var cells = [];

/* ============================================================
   Java-compatible Random (matches java.util.Random exactly)
   ============================================================ */
var JavaRandom = function(seed) {
	const multiplier = 0x5DEECE66Dn;
	const addend = 0xBn;
	const mask = (1n << 48n) - 1n;

	let state = (BigInt(seed) ^ multiplier) & mask;

	function next(bits) {
		state = (state * multiplier + addend) & mask;
		return Number(state >> (48n - BigInt(bits)));
	}

	return {
		nextDouble: function () {
			return ((next(26) << 27) + next(27)) / (1 << 53);
		},
		nextInt: function (bound) {
			return Math.floor(this.nextDouble() * bound);
		}
	};
};

// Default RNG (time-based seed)
var rng = JavaRandom(Date.now()); // default: time-based seed

// Public seed setter
//function setRandomSeed(seed) {
//	rng = JavaRandom(seed);
//	console.log("Random seed set to:", seed);
//}

function setRandomSeed(seed) {
    currentSeed = Number(seed);
    rng = JavaRandom(currentSeed);

    $("#genieSeed").val(currentSeed); // keep UI synced
    console.log("Random seed set to:", currentSeed);
}


;(function($) {

$(document).ready(function() {
	Grid.init();
	NumAllelesPlot.initPlot();
	AlleleFrequencyPlot.initPlot();
});

var colors = [
	"#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999",
	"#714C02", "#01587A", "#024E37", "#726C01", "#003A5E", "#6D2D00", "#752E58", "#4A4A4A",
	"#FFCC94", "#A9DAFF", "#68D7AC", "#FFF360", "#7AB8F9", "#FFA98F", "#FFB1DA", "#CBCBCB"
];

var allele_to_color = function(allele) {
	allele = parseInt(allele);
	if(allele === -2) {
		return "#FFFFFF"
	} else if(allele === -1) {
		return "#000000"
	}
	return colors[allele % colors.length];
}

var Grid = function() {
	//var cells = [];
	var stateCapture = [], simulation, mutationRate = 0;
	var uniqueCells = {}, uniqueMutations = {};
	var numAllelesOverTime = [];
	var numIntervals = 0;
	var alleleFrequencies = {};
	var numBarriers = 0;
	var isRunning = false;

	var theta = 2048 * 0.001; // initial theta

	var update_cell = function(cell_num, allele) {
		if(allele !== undefined) {
			cells[cell_num] = allele;
		}
		var color = allele_to_color(cells[cell_num]);
		$("#cell-" + cell_num).css("background-color", color);
	}

	var getRandomAllele = function() {
		return rng.nextInt(colors.length * 1000);
	}

	var step = function(mutationRate) {
		//Random cell dies
		var cellNum = rng.nextInt(1024);

		//Dead cell is randomly replaced by one of it's neighbors
		var numCellsPerRow = 32, neighbors = [];
		var upper, lower, left, right, upperLeft, upperRight, lowerLeft, lowerRight;
		var randomNeighborIndex, randomNeighbor;

		if (cells[cellNum] == -1) {
			return;
		}

		if (cellNum === 0 || cellNum === 31 || cellNum === 992 || cellNum === 1023) {
			switch(cellNum) {
				case 0:
				// console.log("Corner cell chosen: " + cellNum);
				right = 1;
				lower = 32;
				lowerRight = 33;
				neighbors = [right, lower, lowerRight];
				break;

				case 31:
				// console.log("Corner cell chosen: " + cellNum);
				left = 30;
				lower = 63;
				lowerLeft = 62;
				neighbors = [left, lower, lowerLeft];
				break;

				case 992:
				// console.log("Corner cell chosen: " + cellNum);
				upper = 960;
				upperRight = 961;
				right = 993;
				neighbors = [upper, upperRight, right];
				// alert("Gotcha");
				break;

				case 1023:
				// console.log("Corner cell chosen: " + cellNum);
				upperLeft = 990;
				upper = 991;
				left = 1022;
				neighbors = [upperLeft, upper, left];
				// alert("Gotcha??");
				break;
			}
		}
		else if (cellNum % 32 === 0) {
			// console.log("Left side cell chosen: " + cellNum);
			upper = cellNum - numCellsPerRow;
			upperRight = upper + 1;
			left = cellNum + 1;
			lower = cellNum + numCellsPerRow;
			lowerRight = lower + 1;
			neighbors = [upper, upperRight, left, lower, lowerRight];

		}
		else if (cellNum % 32 === 31) {
			// console.log("Right side cell chosen: " + cellNum);
			upper = cellNum - numCellsPerRow;
			upperLeft = upper - 1;
			left = cellNum - 1;
			lower = cellNum + numCellsPerRow;
			lowerLeft = lower - 1;
			neighbors = [upper, upperLeft, left, lower, lowerLeft];
		}
		else if (cellNum < 32) {
			// console.log("Top cell chosen: " + cellNum);
			left = cellNum - 1;
			right = cellNum + 1;
			lower = cellNum + 32;
			lowerLeft = lower - 1;
			lowerRight = lower + 1;
			neighbors = [left, right, lower, lowerLeft, lowerRight];
		}
		else if (cellNum > 991) {
			// console.log("Bottom cell chosen: " + cellNum);
			left = cellNum - 1;
			right = cellNum + 1;
			upper = cellNum - 32;
			upperLeft = upper - 1;
			upperRight = upper + 1;
			neighbors = [left, right, upper, upperLeft, upperRight];
		}
		else {
			numCellsPerRow = 32;
			upper = cellNum - numCellsPerRow;
			lower = cellNum + numCellsPerRow;
			left = cellNum - 1;
			right = cellNum + 1;
			upperLeft = upper - 1;
			upperRight = upper + 1;
			lowerLeft = lower - 1;
			lowerRight = lower + 1;
			neighbors = [upperLeft, upper, upperRight, left, right, lowerLeft, lower, lowerRight];
		}

		//New cell will mutate with probability given by mutation rate
		var rand = rng.nextDouble();
		if (rand < mutationRate) {
			update_cell(cellNum, getRandomAllele());
		}
		else {
			for(let i = 0; i < 20; i ++) {
				randomNeighborIndex = rng.nextInt(neighbors.length);
				randomNeighbor = neighbors[randomNeighborIndex];
				// make sure that it's not a dead cell
				if(cells[randomNeighbor] >= 0) {
					update_cell(cellNum, cells[randomNeighbor]);
					break;
				}
			}
		}
	};


	var drawGrid = function() {
		// count the number of occurrences of each unique element of cells
		// https://stackoverflow.com/a/43673826
		var counts = cells.reduce( (acc, o) => (acc.set(o, (acc.get(o) || 0) + 1), acc), new Map() );

		// remove empty and barrier cells
		counts.delete(-1);
		counts.delete(-2);

		var num_alleles = counts.size;
		numAllelesOverTime.push([numIntervals, num_alleles]);

		var num_alive = Array.from(counts.values()).reduce( (acc, o) => (acc + o), 0);

		for (const [key, value] of counts.entries()) {

			if (typeof(alleleFrequencies[key]) == "undefined") {
				alleleFrequencies[key.toString()] = {color:allele_to_color(key), data:[]}
			}
			alleleFrequencies[key.toString()].data.push([numIntervals, value/num_alive]);
		}

		var lines = Object.keys(alleleFrequencies).map((k) => (alleleFrequencies[k]));

		NumAllelesPlot.update(numAllelesOverTime);
		AlleleFrequencyPlot.update(lines, numIntervals);
	};

	var generateStatistics = function() {
		var uniqueCells = {};
		var uniqueMutations = {};
		var numAlleles = 0;
		var numActiveMutations = 0;
		for (var i = 0; i < cells.length; i++) {
			uniqueCells[cells[i]] = true;
			if (cells[i] > 1024) {
				uniqueMutations[cells[i]] = true;
			}
		}
		for (var i in uniqueCells) {
			numAlleles++;
		}
		for (var i in uniqueMutations) {
			numActiveMutations++;
		}
		$("#numAlleles").html("Number of alleles: " + numAlleles);
	};

	var handleSeedButton = function () {
    	$("#genieSetSeed").button().click(function () {
        	const seedVal = $("#genieSeed").val();

        	if (seedVal === "" || !Number.isInteger(Number(seedVal))) {
            	alert("Please enter a valid integer seed");
            	return;
        	}

        	setRandomSeed(Number(seedVal));
		});
	};

	var handleStartButton = function() {
		$("#genieStart").button({label: "Start", icon: "fa fa-play"}).click(function() {
			if (!isRunning) {
				mutationRate = $("#genieMutationRate").val();
				if (mutationRate === null || mutationRate === "") {
					alert("Please enter a mutation rate");
				}
				else if (isNaN(mutationRate)) {
					alert("Mutation rate must be numeric");
				}
				else if (mutationRate > 1 || mutationRate < 0) {
					alert("Mutation rate must be between 0 and 1");
				}
				else {
					simulation = setInterval(function() {
						//$("#exportButton").addClass("disabled");
						$("#genieStart").button("option", "label", "Pause");
						$("#genieStart").button("option", "icon", "fa fa-pause");						
						$("#genieReset").button("option", "disabled", true);
						runSimulation();
						drawGrid();
						numIntervals++;
						isRunning = true;
					}, 200);
				}
			}
			else {
				$("#genieStart").button("option", "label", "Start");
				$("#genieStart").button("option", "icon", "fa fa-play");
				$("#genieReset").button("option", "disabled", false);
				//$("#exportButton").removeClass("disabled");
				//createExportLink();

				clearInterval(simulation);
				isRunning = false;
				generateStatistics();
			}
		});
	};

	var runSimulation = function() {
		for (var i = 0; i < 2000; i++) {
			step(mutationRate);
		}
	};

	var handleResetButton = function() {
		$("#genieReset").button({label: "Reset", icon: "fa fa-refresh"}).click(function() {
			cells = [];

			$(".show-cell").each(function(index) {
				if (index !== 0 && rng.nextDouble() < index/(index + theta)) {
					var cellNum = rng.nextInt((index - 1));
					while (cells[cellNum] == -1) {
						cellNum = rng.nextInt((index - 1));
					}
					cells.push(cells[cellNum]);
				}
				else {
					var allele = getRandomAllele();
					cells.push(allele);
				}
				$(this).attr("id", "cell-" + index);
				update_cell(index);
			});
			generateStatistics();
			alleleFrequencies = [];
			numIntervals = 0;
			numBarriers = 0;
			numAllelesOverTime = [];
			NumAllelesPlot.initPlot();
			AlleleFrequencyPlot.initPlot();
			//createExportLink();
		});
	};

	var createExportLink = function() {
		$("#exportButton").unbind("click");
		$("#exportButton").click(function() {
			var returnObj = {
				cells: cells
			};
			var data = "text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(returnObj));
			window.open("data: " + data);
		});
		//$("#exportButton").html("<a id='exportLink' href='data:" + data + "' download='data.json'>Export Grid</a>");
	};

	var handleBarrier = function() {
		$(".show-cell").mousedown(function(event) {
			if (!event.shiftKey) {
				event.preventDefault();
				var cellNum = parseInt($(this).attr("id").slice(5));
				if (cells[cellNum] !== -1) {
					cells[cellNum] = -1;
					update_cell(cellNum);
					numBarriers++;
				}
				else {
					cells[cellNum] = -2;
					update_cell(cellNum);
					numBarriers--;
				}
				$(".show-cell").mouseover(function(e) {
					if (!e.shiftKey) {
						var cellNum = parseInt($(this).attr("id").slice(5));
						if (cells[cellNum] !== -1) {
							cells[cellNum] = -1;
							update_cell(cellNum);
							numBarriers++;
						}
						else {
							cells[cellNum] = -2;
							update_cell(cellNum);
							numBarriers--;
						}
					}
				});
			}
		});
		$(document).mouseup(function() {
			$(".show-cell").unbind("mouseover");
			//createExportLink();
		});
	};

	var handleBarrierTemplate = function() {
		$("#barrierTemplateButton").click(function() {
			for (var i = 0; i < 1024; i++) {
				if (barrierGrid.cells[i] == -1) {
					cells[i] = -1;
					update_cell(i);
					numBarriers++;
				}
			}
		});
	};

	var handleForcedMutation = function() {
		$(".show-cell").mousedown(function(event) {
			if (event.shiftKey) {
				event.preventDefault();
				var cellNum = parseInt($(this).attr("id").slice(5));
				console.log("New mutation at cell %d", cellNum);

				var newAllele = getRandomAllele();
				cells[cellNum] = newAllele;
				update_cell(cellNum);
				
				$(".show-cell").mouseover(function(event) {
					if (event.shiftKey) {
						var cellNum = parseInt($(this).attr("id").slice(5));
						cells[cellNum] = newAllele;
						update_cell(cellNum);
					}
				});
			}
		});
		$(document).mouseup(function() {
			$(".show-cell").unbind("mouseover");
			//createExportLink();
		});
	};

	var handleHelpButton = function() {
		$("#genieHelp").button().click(function() {
			window.open("https://github.com/benhroos/thesis/blob/master/README.md");
		});
	};

	var resizePlot = function() {
		$("#numAllelesPlot").height($("#genieGrid").height()/2.25);
		$("#alleleFrequencyPlot").height($("#genieGrid").height()/2.25);
		$(window).resize(function() {
			$("#numAllelesPlot").height($("#genieGrid").height()/2.25);
			$("#alleleFrequencyPlot").height($("#genieGrid").height()/2.25);
		});
	};

	var init = function() {
		var gridRowHTML = $("#genieGrid").html();
		for (var count = 0; count < 31; count++) {
			$("#genieGrid").append(gridRowHTML);
		}

		$(".show-cell").each(function(index) {
			if (index !== 0 && rng.nextDouble() < index/(index + theta)) {
				var cellNum = rng.nextInt((index - 1));
				while (cells[cellNum] == -1) {
					cellNum = rng.nextInt((index - 1));
				}
				cells.push(cells[cellNum]);
			}
			else {
				var allele = getRandomAllele();
				cells.push(allele);
			}
			$(this).attr("id", "cell-" + index);
			update_cell(index);
		});
		//createExportLink();
		handleStartButton();
		handleResetButton();
		handleSeedButton();
		handleBarrier();
		handleBarrierTemplate();
		handleForcedMutation();
		handleHelpButton();
		resizePlot();
		generateStatistics();
	};

	return {
		init: init,
		getCells: function() {
			return cells;
		},
		getCapture: function() {
			return stateCapture;
		},
		stats: generateStatistics,
		getUniqueAlleles: function() {
			return uniqueCells;
		},
		getAlleleFrequencies: function() {
			return alleleFrequencies;
		},
		getNumIntervals: function() {
			return numIntervals;
		}
	}
}();

var NumAllelesPlot = (function() {
    var plot;
    function init() {
        plot = $.plot("#numAllelesPlot", [[]], {
            axisLabels: {
                show: true,
            },
            xaxis: {
                axisLabel: "Generation",
                min: 0,
                max: 100,
                tickSize: 10
            },
            yaxis: {
                axisLabel: "Alleles",
                min: 0,
                max: 16,
                tickSize: 2
            }
        });
    }

    function updatePlot(numAllelesOverTime) {
        var gridLength = numAllelesOverTime.length;
        var xMin, xMax, maxPan;
        if (gridLength > 100) {
            xMin = gridLength - 100;
            xMax = gridLength - 1;
            maxPan = gridLength - 1;
        }
        else {
            xMin = 0;
            xMax = 100;
            maxPan = 100;
        }

        plot = $.plot("#numAllelesPlot", [numAllelesOverTime], {
            axisLabels: {
                show: true,
            },
            xaxis: {
                axisLabel: "Generation",
                min: xMin,
                max: xMax,
                tickSize: 10,
                panRange: [0, maxPan]
            },
            yaxis: {
                axisLabel: "Alleles",
                min: 0,
                max: 16,
                tickSize: 2,
                panRange: [0, 16]
            },
            pan: {
                interactive: true
            }
        });

        // plot.setupGrid();
        // plot.draw();
    }

    return {
        initPlot: init,
        update: updatePlot,
        getPlot: function() {
            return plot;
        },
        redrawPlot: function() {
            plot.draw();
        }
    }
})();

var AlleleFrequencyPlot = (function() {
    var plot;
    function init() {
        plot = $.plot("#alleleFrequencyPlot", [[]], {
            axisLabels: {
                show: true,
            },
            xaxis: {
                axisLabel: "Generation",
                min: 0,
                max: 100,
                tickSize: 10
            },
            yaxis: {
                axisLabel: "Allele Frequency",
                min: 0,
                max: 1,
                tickSize: .1
            }
        });
    }

    function updatePlot(alleleFrequencies, numIntervals) {
        var gridLength = numIntervals;
        var xMin, xMax, maxPan;
        if (gridLength > 100) {
            xMin = gridLength - 100;
            xMax = gridLength - 1;
            maxPan = gridLength - 1;
        }
        else {
            xMin = 0;
            xMax = 100;
            maxPan = 100;
        }

        plot = $.plot("#alleleFrequencyPlot", alleleFrequencies, {
            axisLabels: {
                show: true,
            },
            xaxis: {
                axisLabel: "Generation",
                min: xMin,
                max: xMax,
                tickSize: 10,
                panRange: [0, maxPan]
            },
            yaxis: {
                axisLabel: "Allele Frequency",
                min: 0,
                max: 1,
                tickSize: .1,
                panRange: [0, 1]
            },
            pan: {
                interactive: true
            }
        });
    }

    return {
        initPlot: init,
        update: updatePlot,
        getPlot: function() {
            return plot;
        }
    }
})();

var barrierGrid = {"cells":[{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#cab2d6","allele":49,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#cab2d6","allele":49,"mutationNumber":-1},{"color":"#cab2d6","allele":49,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#cab2d6","allele":49,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#cab2d6","allele":49,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#cab2d6","allele":49,"mutationNumber":-1},{"color":"#cab2d6","allele":49,"mutationNumber":-1},{"color":"#cab2d6","allele":49,"mutationNumber":-1},{"color":"#cab2d6","allele":49,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#cab2d6","allele":49,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#e31a1c","allele":8,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#1f78b4","allele":1,"mutationNumber":-1},{"color":"#ff7f00","allele":33,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#b2df8a","allele":7,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#000000","allele":-1,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1},{"color":"#33a02c","allele":3,"mutationNumber":-1}],"numAlleles":7,"numMutations":0};

})(jQuery);
