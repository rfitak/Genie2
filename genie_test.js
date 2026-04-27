// Copyright 2015 Ben H. Roos <beenhroos@gmail.com>
// Copyright 2016 Reed A. Cartwright <reed@cartwrig.ht>

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
var rng = JavaRandom(Date.now());

// Public seed setter
function setRandomSeed(seed) {
    rng = JavaRandom(seed);
    console.log("Random seed set to:", seed);
}

;(function($) {
$(document).ready(function() {
    Grid.init();
    NumAllelesPlot.initPlot();
    AlleleFrequencyPlot.initPlot();
});

var colors = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999",
    "#714C02", "#01587A", "#024E37", "#726C01", "#003A5E", "#6D2D00", "#752E58", "#4A4A4A"
];

var allele_to_color = function(allele) {
    allele = parseInt(allele);
    if (allele === -2) return "#FFFFFF";
    if (allele === -1) return "#000000";
    return colors[allele % colors.length];
};

var Grid = function() {
    var stateCapture = [], simulation, mutationRate = 0;
    var uniqueCells = {}, alleleFrequencies = {};
    var numAllelesOverTime = [], numIntervals = 0;
    var numBarriers = 0, isRunning = false;
    var theta = 2048 * 0.001;

    var update_cell = function(cell_num, allele) {
        if (allele !== undefined) cells[cell_num] = allele;
        $("#cell-" + cell_num)
            .css("background-color", allele_to_color(cells[cell_num]));
    };

    var getRandomAllele = function() {
        return rng.nextInt(colors.length * 1000);
    };

    var step = function(mutationRate) {
        var cellNum = rng.nextInt(1024);
        if (cells[cellNum] === -1) return;

        var neighbors = [];
        var row = Math.floor(cellNum / 32);
        var col = cellNum % 32;

        for (let dr = -1; dr <= 1; dr++) {
            for (let dc = -1; dc <= 1; dc++) {
                if (dr === 0 && dc === 0) continue;
                let nr = row + dr, nc = col + dc;
                if (nr >= 0 && nr < 32 && nc >= 0 && nc < 32) {
                    neighbors.push(nr * 32 + nc);
                }
            }
        }

        if (rng.nextDouble() < mutationRate) {
            update_cell(cellNum, getRandomAllele());
        } else {
            for (let i = 0; i < 20; i++) {
                let n = neighbors[rng.nextInt(neighbors.length)];
                if (cells[n] >= 0) {
                    update_cell(cellNum, cells[n]);
                    break;
                }
            }
        }
    };

    var drawGrid = function() {
        var counts = cells.reduce(
            (acc, o) => (acc.set(o, (acc.get(o) || 0) + 1), acc),
            new Map()
        );
        counts.delete(-1);
        counts.delete(-2);

        numAllelesOverTime.push([numIntervals, counts.size]);
        var alive = Array.from(counts.values()).reduce((a, b) => a + b, 0);

        counts.forEach((v, k) => {
            if (!alleleFrequencies[k]) {
                alleleFrequencies[k] = {
                    color: allele_to_color(k),
                    data: []
                };
            }
            alleleFrequencies[k].data.push([numIntervals, v / alive]);
        });

        NumAllelesPlot.update(numAllelesOverTime);
        AlleleFrequencyPlot.update(
            Object.values(alleleFrequencies),
            numIntervals
        );
    };

    var runSimulation = function() {
        for (var i = 0; i < 2000; i++) step(mutationRate);
    };

    var handleStartButton = function() {
        $("#genieStart").button().click(function() {
            if (!isRunning) {
                mutationRate = parseFloat($("#genieMutationRate").val());
                simulation = setInterval(function() {
                    runSimulation();
                    drawGrid();
                    numIntervals++;
                    isRunning = true;
                }, 200);
            } else {
                clearInterval(simulation);
                isRunning = false;
            }
        });
    };

    var handleResetButton = function() {
        $("#genieReset").button().click(function() {
            cells = [];
            $(".show-cell").each(function(index) {
                if (index !== 0 && rng.nextDouble() < index / (index + theta)) {
                    let cellNum = rng.nextInt(index);
                    cells.push(cells[cellNum]);
                } else {
                    cells.push(getRandomAllele());
                }
                $(this).attr("id", "cell-" + index);
                update_cell(index);
            });
            alleleFrequencies = {};
            numIntervals = 0;
            numAllelesOverTime = [];
            NumAllelesPlot.initPlot();
            AlleleFrequencyPlot.initPlot();
        });
    };

    var init = function() {
        var rowHTML = $("#genieGrid").html();
        for (var i = 0; i < 31; i++) $("#genieGrid").append(rowHTML);

        $(".show-cell").each(function(index) {
            let allele = (index !== 0 && rng.nextDouble() < index / (index + theta))
                ? cells[rng.nextInt(index)]
                : getRandomAllele();
            cells.push(allele);
            $(this).attr("id", "cell-" + index);
            update_cell(index);
        });

        handleStartButton();
        handleResetButton();
    };

    return { init: init };
}();

/* === Plot code unchanged === */
})(jQuery);
