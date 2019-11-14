#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create samplot vcf commands to execute and generate
companion HTML image browser.

Note: additional arguments are passed through to samplot plot
"""
from __future__ import print_function

import json
import operator
import os
import random
import sys

import pysam

try:
    from shlex import quote
except ImportError:
    from pipes import quote


HERE = os.path.dirname(__file__)
HTML = """<!DOCTYPE html>
<html lang='en'>
<head>
    <meta charset='utf-8'>
    <title>samplot</title>

    <script src="https://cdnjs.cloudflare.com/ajax/libs/d3/5.9.2/d3.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/crossfilter2/1.4.7/crossfilter.min.js" type='text/javascript'></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/dc/3.0.12/dc.min.js" type='text/javascript'></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.4.1/jquery.min.js" type='text/javascript'></script>
    <!-- 1.10.16 due to export fixes in latest -->
    <script src="https://cdn.datatables.net/1.10.16/js/jquery.dataTables.js" type='text/javascript'></script>
    <script src="https://cdn.datatables.net/1.10.16/js/dataTables.bootstrap4.min.js" type='text/javascript'></script>
    <!-- export buttons -->
    <script src="https://cdn.datatables.net/buttons/1.5.6/js/dataTables.buttons.min.js" type='text/javascript'></script>
    <script src="https://cdn.datatables.net/buttons/1.5.6/js/buttons.html5.min.js" type='text/javascript'></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/4.3.1/js/bootstrap.bundle.min.js" type='text/javascript'></script>

    <link href="https://cdn.datatables.net/buttons/1.5.6/css/buttons.dataTables.min.css" rel="stylesheet" type="text/css"/>
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.8.2/css/all.min.css" rel='stylesheet' type='text/css'>
    <link href="https://cdn.datatables.net/v/bs4/dt-1.10.18/sl-1.3.0/datatables.min.css" rel="stylesheet" type="text/css"/>
    <link href="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/4.3.1/css/bootstrap.min.css" rel='stylesheet' type='text/css'>
    <link href="https://cdnjs.cloudflare.com/ajax/libs/dc/3.0.12/dc.min.css" rel='stylesheet' type='text/css'>

    <style type="text/css">
        #filter-menu .dropdown-menu { min-height: 100px; max-height: 100vh; overflow-y: auto; overflow-x: hidden; background-color: #edf0f2; }
        span.no-show { display: none; }
        span.show-ellipsis:after { content: "..."; }
        .datatable-info { font-size: .9em; }
        #variant-table_info { padding-top: 8px; }
        table.dataTable thead th.sorting:after,
        table.dataTable thead th.sorting_asc:after,
        table.dataTable thead th.sorting_desc:after,
        table.dataTable thead th.sorting:before,
        table.dataTable thead th.sorting_asc:before,
        table.dataTable thead th.sorting_desc:before {
            font-family: FontAwesome !important;
        }
        .modal-content { width: 610px; }
        h7 { font-size: .95rem; }
    </style>
</head>

<body>
    <nav class="navbar navbar-dark bg-dark p-0 pl-2">
        <a class="navbar-brand text-light p-0" href="https://github.com/ryanlayer/samplot">samplot</a>
    </nav>

    <div class="modal fade" id="filter-modal" tabindex="-1" role="dialog" aria-labelledby="filter-modal" aria-hidden="true">
        <div class="modal-dialog" role="document">
            <div class="modal-content">
                <div class="modal-header">
                    <div class="flex-column">
                        <h5 class="modal-title" id="filter-modal">Filters</h5>
                        <h7 class="pl-2 text-secondary" id="variant-count">
                            <a href="javascript:dc.filterAll(); dc.renderAll();">Reset All</a>
                        </h7>
                    </div>
                </div>
                <div class="modal-body">
                    <div class="container">
                        <div class="row pt-2">
                            <div class="col">
                                <h5>Sample</h5>
                            </div>
                        </div>
                        <div class="row pb-3">
                            <div class="col-12">
                                <div id="sample-search"></div>
                            </div>
                        </div>
                        <div class="row" id="nsamples-chart">
                            <div class="col-4">
                                <h5># of Samples</h5>
                            </div>
                            <div class="col-8 text-right">
                                <span class="reset text-muted" style="display: none;">[<span class="filter"></span>]</span>
                                <a class="reset" href="javascript:nsamplesChart.filterAll();dc.redrawAll();" style="display: none;">Reset</a>
                            </div>
                        </div>
                        <div class="row" id="size-chart">
                            <div class="col-4">
                                <h5>Size</h5>
                            </div>
                            <div class="col-8 text-right">
                                <span class="reset text-muted" style="display: none;">[<span class="filter"></span>]</span>
                                <a class="reset" href="javascript:sizeChart.filterAll();dc.redrawAll();" style="display: none;">Reset</a>
                            </div>
                        </div>
                        <div class="row" id="type-chart">
                            <div class="col-4">
                                <h5>SV Type</h5>
                            </div>
                            <div class="col-8 text-right">
                                <span class="reset text-muted" style="display: none;">[<span class="filter"></span>]</span>
                                <a class="reset" href="javascript:typeChart.filterAll();dc.redrawAll();" style="display: none;">Reset</a>
                            </div>
                        </div>
                        <div class="row" id="chrom-chart">
                            <div class="col-4">
                                <h5>Chromosome</h5>
                            </div>
                            <div class="col-8 text-right">
                                <span class="reset text-muted" style="display: none;">[<span class="filter"></span>]</span>
                                <a class="reset" href="javascript:chromChart.filterAll();dc.redrawAll();" style="display: none;">Reset</a>
                            </div>
                        </div>
                        <div class="row" id="overlaps-chart" hidden>
                            <div class="col-4">
                                <h5>SV Overlaps</h5>
                            </div>
                            <div class="col-8 text-right">
                                <span class="reset text-muted" style="display: none;">[<span class="filter"></span>]</span>
                                <a class="reset" href="javascript:overlapsChart.filterAll();dc.redrawAll();" style="display: none;">Reset</a>
                            </div>
                        </div>
                    </div>
                </div>
                <div class="modal-footer">
                    <button type="button" class="btn btn-outline-secondary" data-dismiss="modal" onclick="javascript:dc.filterAll(); dc.renderAll();" title="Clear selection and close">Cancel</button>
                    <button type="button" class="btn btn-primary" data-dismiss="modal" title="Apply filters">Apply</button>
                </div>
            </div>
        </div>
    </div>

    <div class="container-fluid bg-light">
        <div class="row" id="variant-table-placeholder" style="height:435px">
            <div class="col-12">
                <div class="d-flex justify-content-center align-items-center text-muted h-100">
                    <div class="d-flex flex-column">
                        <i class="fas fa-7x fa-table"></i>
                    </div>
                </div>
            </div>
        </div>

        <div class="row pb-1" id="variant-table-div" hidden>
            <div class="col-12">
                <div class="table-responsive">
                    <table id="variant-table" class="table table-hover table-sm table-bordered display nowrap" width="100%"></table>
                </div>
            </div>
        </div>
    </div>

    <div class="container-fluid">
        <div class="row" id="samplot-area-placeholder">
            <div class="col-12">
                <div style="height:415px">
                    <div class="d-flex justify-content-center align-items-center text-muted h-100">
                        <div class="d-flex flex-column">
                            <i class="fas fa-7x fa-chart-bar"></i>
                            <ul class="small p-0 m-3">
                                <li>
                                    Select a variant </br>from the table
                                </li>
                            </ul>
                        </div>
                    </div>
                </div>
            </div>
        </div>
        <div class="row" id="samplot-area" hidden>
            <div class="col-12">
                <img id="samplot-img" width="100%" src=""/>
            </div>
        </div>
    </div>
</body>

<script>
const data = [DATA]
const plot_type = ".[PLOT_TYPE]"
const annotation = [GFF]
const denovo = [DENOVO]

dc.config.defaultColors(d3.schemeSet1)

// plot constraints
const plotw = 585
const ploth = 150

// table filters
var searchInput = dc.textFilterWidget("#sample-search")
var nsamplesChart = dc.barChart("#nsamples-chart")
var sizeChart = dc.barChart("#size-chart")
var typeChart = dc.barChart("#type-chart")
var chromChart = dc.barChart("#chrom-chart")
var overlapsChart
// shows filter impact in modal header
var variantCount = dc.dataCount("#variant-count")

// used to access filtered table data
var chromDimension
// datatables obj
var variant_table
// crossfilter obj
var ndx

$('#filter-modal').on('hidden.bs.modal', function () {
    update_table()
})

const table_click = (row) => {
    d3.select('#samplot-area-placeholder')
        .property("hidden", true)
    d3.select('#samplot-area')
        .property("hidden", false)

    // load the image
    var img = row.svtype + "_" + row.chrom + "_" + row.start + "_" + row.end + plot_type
    d3.select("#samplot-img")
        .property("src", img)
}

function add_download_actions() {
    $('#copy-button-download').on('click', function() {
        variant_table.button('.buttons-copy').trigger()
    })
    $('#csv-button-download').on('click', function() {
        variant_table.button('.buttons-csv').trigger()
    })
}

function build_table(data) {
    // hide the placeholder and show the datatable
    d3.select('#variant-table-placeholder').property("hidden", true)
    d3.select('#variant-table-div').property("hidden", false)

    let cols = [
        {data: 'chrom', title: 'Chrom'},
        {data: 'start', title: 'Start'},
        {data: 'end', title: 'End'},
        {data: 'svlength', title: 'Size'},
        {data: 'svtype', title: 'SV Type'},
        {data: 'nsamples', title: '# of Samples'},
        {data: 'samples', title: 'Samples'},
    ]
    if (annotation) {
        d3.select('#overlaps-chart').property("hidden", false)
        cols.push({data: 'overlaps', title: 'Overlaps'})
    }
    if (denovo) {
        cols.push({data: 'dn', title: 'De novo'})
    }

    variant_table = $("#variant-table").DataTable({
        data: data,
        columns: cols,
        scrollX: true,
        scrollCollapse: true,
        paging: true,
        pagingType: "full",
        pageLength: 10,
        info: true,
        buttons: [
            'copyHtml5', 'csvHtml5'
        ],
        fnInfoCallback: (oSettings, iStart, iEnd, iMax, iTotal, sPre) => {
            return `
            <span class="datatable-info">
                <span class="pr-2">Showing <b>${iStart}</b> - <b>${iEnd}</b> of <b>${iTotal}</b> records</span>
                <button type="button" class="btn btn-primary btn-sm" data-toggle="modal" data-target="#filter-modal" title="Show filters">
                    <span class="fas fa-filter"></span>
                </button>
                <span class="dropup">
                    <button type="button" class="btn btn-sm btn-primary dropdown-toggle" id="download-menu" title="Save table" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
                        <span class="fas fa-save"></span>
                    </button>
                    <span class="dropdown-menu" aria-labelledby="download-menu">
                        <h6 class="dropdown-header">Save ${iTotal} rows as:</h6>
                        <button class="dropdown-item" type="button" id="csv-button-download">
                            CSV
                        </button>
                        <button class="dropdown-item" type="button" id="copy-button-download">
                            Copy
                        </button>
                    </span>
                </span>
            </span>
            `
        },
        columnDefs: [
            {
                targets: (annotation ? [0,1,2,3,4,5,7] : [0,1,2,3,4,5]),
                width: '15%'
            },
            // https://datatables.net/blog/2016-02-26
            {
                targets: 6,
                render: function(data, type, row) {
                    if (type === 'display' && data != null) {
                        data = data.replace(/<(?:.|\\n)*?>/gm, '');
                        if(data.length > 40) {
                            return '<span class=\"show-ellipsis\">' + data.substr(0, 40) + '</span><span class=\"no-show\">' + data.substr(40) + '</span>';
                        } else {
                            return data;
                        }
                    } else {
                        return data;
                    }
                }
            }
        ],
        // search is applied using crossfilter
        searching: false,
        lengthChange: false,
        order: [[0, 'asc'], [1, 'asc']]
    })

    // register table clicks on sample_column
    variant_table.on('click', 'tr', function () {
        if ( $(this).hasClass('selected') ) {
            // do nothing for now
        }
        else {
            // de-select all
            variant_table.$('tr.selected').removeClass('selected')
            // add select to clicked row
            $(this).addClass('selected')
            // update samplot area
            table_click(variant_table.rows('.selected').data()[0])
        }
    })

    add_download_actions()
}

function update_table() {
    variant_table.clear()
    variant_table.rows.add(chromDimension.top(Infinity))
    variant_table.draw()
    // gives function to new buttons each time table is redrawn
    add_download_actions()
}

function remove_empty_bins(source_group) {
    return {
        all:function () {
            return source_group.all().filter(function(d) {
                return d.value != 0
            })
        }
    }
}

// https://jsfiddle.net/gordonwoodhull/g34Ldwaz/8/
// https://github.com/dc-js/dc.js/issues/348
function index_group(group) {
    return {
        all: function() {
            return group.all().map(function(kv, i) {
                return {key: i, value: kv.value}
            })
        }
    }
}

$(document).ready(function() {

    ndx = crossfilter(data)
    var all = ndx.groupAll()

    chromDimension = ndx.dimension((d) => { return d.chrom })
    build_table(chromDimension.top(Infinity))
    var chromGroup = chromDimension.group().reduceCount()
    var nonEmptyChromGroup = remove_empty_bins(chromGroup)

    var searchDimension = ndx.dimension(function(d) {
        return d.samples
    })
    searchInput
        .dimension(searchDimension)
        .on('renderlet', function() {
            d3.selectAll(".dc-text-filter-input")
                .classed("form-control", true)
            d3.selectAll("#sample-search.dc-chart")
                .classed("col-12", true)
        })

    var sizeDimension = ndx.dimension(function(d) {
        var round
        if (d.svlength < 100) {
            round = 100
        } else if (d.svlength < 1000) {
            round = 100
        } else if (d.svlength < 10000) {
            round = 1000
        } else if (d.svlength < 100000) {
            round = 10000
        } else if (d.svlength < 1000000) {
            round = 100000
        } else if (d.svlength < 10000000) {
            round = 1000000
        } else {
            round = 10000000
        }
        return Math.round(d.svlength / round) * round
    })
    var sizeGroup = sizeDimension.group().reduceCount()
    var nonEmptySizeGroup = remove_empty_bins(sizeGroup)
    // for brushing, need to track keys at numeric indexes
    var sizeKeys = nonEmptySizeGroup.all().map(dc.pluck('key')).slice()

    var typeDimension = ndx.dimension((d) => { return d.svtype })
    var typeGroup = typeDimension.group().reduceCount()
    var nonEmptyTypeGroup = remove_empty_bins(typeGroup)

    var nsamplesDimension = ndx.dimension((d) => { return d.nsamples })
    var nsamplesDimension = ndx.dimension(function(d) {
        var round
        if (d.nsamples < 10) {
            round = 1
        } else if (d.nsamples < 100) {
            round = 10
        } else if (d.nsamples < 1000) {
            round = 100
        } else if (d.nsamples < 10000) {
            round = 1000
        } else {
            round = 10000
        }
        return Math.round(d.nsamples / round) * round
    })
    var nsamplesGroup = nsamplesDimension.group().reduceCount()
    var nonEmptyNsamplesGroup = remove_empty_bins(nsamplesGroup)
    var nsamplesKeys = nonEmptyNsamplesGroup.all().map(dc.pluck('key')).slice()

    // number of samples
    nsamplesChart
        .width(plotw).height(ploth).gap(1)
        .margins({top: 10, right: 50, bottom: 30, left: 40})
        .x(d3.scaleLinear().domain([0, nsamplesKeys.length]))
        .round(Math.floor)
        .brushOn(true)
        .elasticX(true)
        .dimension(nsamplesDimension)
        .group(index_group(nonEmptyNsamplesGroup))
        .elasticY(true)
        .yAxisLabel('Count')
        .filterPrinter(function (filters) {
            var filter = filters[0]
            return nsamplesKeys[filter[0]] + ' - ' + nsamplesKeys[filter[1]]
        })
    // limit the number of labels along x-axis
    nsamplesChart.xAxis().ticks(20)
    nsamplesChart.yAxis().ticks(5)
    // update labels from keys
    nsamplesChart.xAxis().tickFormat(function(v) {
        return nsamplesKeys[v]
    })
    nsamplesChart.filterHandler(function (dimension, filters) {
        if (filters.length === 0) {
            // the empty case (no filtering)
            dimension.filter(null)
        } else {
            dimension.filterRange([nsamplesKeys[filters[0][0]], nsamplesKeys[filters[0][1]]])
        }
        return filters
    })

    // SV length
    sizeChart
        .width(plotw).height(ploth).gap(1)
        .margins({top: 10, right: 50, bottom: 30, left: 40})
        .x(d3.scaleLinear().domain([0, sizeKeys.length]))
        .round(Math.floor)
        .brushOn(true)
        .elasticX(true)
        .dimension(sizeDimension)
        .group(index_group(nonEmptySizeGroup))
        .elasticY(true)
        .yAxisLabel('Count')
        .filterPrinter(function (filters) {
            var filter = filters[0]
            return sizeKeys[filter[0]] + ' - ' + sizeKeys[filter[1]]
        })
        // adds left padding to plots inside filtering panel
        .on('renderlet', function() {
            d3.selectAll("svg")
                .classed("pl-3", true)
        })
    // limit the number of labels along x-axis
    sizeChart.xAxis().ticks(10)
    sizeChart.yAxis().ticks(5)
    // update labels from keys
    sizeChart.xAxis().tickFormat(function(v) {
        return sizeKeys[v]
    })
    // update the status format for this chart
    sizeChart.filterHandler(function (dimension, filters) {
        if (filters.length === 0) {
            // the empty case (no filtering)
            dimension.filter(null)
        } else {
            dimension.filterRange([sizeKeys[filters[0][0]], sizeKeys[filters[0][1]]])
        }
        return filters
    })

    // sv type
    typeChart
        .width(plotw).height(ploth).gap(1)
        .margins({top: 10, right: 50, bottom: 30, left: 40})
        .x(d3.scaleBand())
        .xUnits(dc.units.ordinal)
        .elasticX(true)
        .elasticY(true)
        .dimension(typeDimension)
        .group(nonEmptyTypeGroup)
        .yAxisLabel('Count')
    typeChart.yAxis().ticks(5)

    // chromosome
    chromChart
        .width(plotw).height(ploth).gap(1)
        .margins({top: 10, right: 50, bottom: 30, left: 40})
        .x(d3.scaleBand())
        .xUnits(dc.units.ordinal)
        .yAxisLabel('Count')
        .elasticX(true)
        .elasticY(true)
        .dimension(chromDimension)
        .group(nonEmptyChromGroup)
        .ordering((d) => {
            v = parseInt(d.key)
            if (v) {
                return v
            } else {
                return d.key
            }
        })
    chromChart.yAxis().ticks(5)

    // overlaps
    if (annotation) {
        var overlapsDimension = ndx.dimension((d) => { return d.overlaps })
        var overlapsGroup = overlapsDimension.group().reduceCount()
        var nonEmptyOverlapsGroup = remove_empty_bins(overlapsGroup)
        overlapsChart = dc.barChart("#overlaps-chart")
        overlapsChart
            .width(plotw).height(ploth).gap(1)
            .margins({top: 10, right: 50, bottom: 30, left: 40})
            .x(d3.scaleBand())
            .xUnits(dc.units.ordinal)
            .elasticX(true)
            .elasticY(true)
            .dimension(overlapsDimension)
            .group(nonEmptyOverlapsGroup)
            .yAxisLabel('Count')
        overlapsChart.yAxis().ticks(5)
    }

    variantCount
        .crossfilter(ndx)
        .groupAll(all)
        // (_optional_) `.html` sets different html when some records or all records are selected.
        // `.html` replaces everything in the anchor with the html given using the following function.
        // `%filter-count` and `%total-count` are replaced with the values obtained.
        .html({
            some: '<strong>%filter-count</strong> selected out of <strong>%total-count</strong> records' +
                ' | <a href=\\'javascript:dc.filterAll(); dc.renderAll();\\'>Reset All</a>',
            all: '<strong>%total-count</strong> records'
        });

    dc.renderAll()
})

</script>
</html>
"""
cmp_lookup = {
    ">": operator.gt,  # e.g. DHFC < 0.5
    "<": operator.lt,
    "<=": operator.le,
    ">=": operator.ge,
    "==": operator.eq,
    "contains": operator.contains,  # e.g. CSQ contains HIGH
    "exists": lambda a, b: True,  # e.g. exists smoove_gene
}


class Sample(object):
    __slots__ = [
        "family_id",
        "id",
        "paternal_id",
        "maternal_id",
        "mom",
        "dad",
        "kids",
        "i",
    ]

    def __init__(self, line):
        toks = line.rstrip().split()
        self.family_id = toks[0]
        self.id = toks[1]
        self.paternal_id = toks[2]
        self.maternal_id = toks[3]
        self.kids = []
        self.i = -1  # index in the vcf.

    def __repr__(self):
        return "Sample(id:{id},paternal_id:{pid},maternal_id:{mid})".format(
            id=self.id, pid=self.paternal_id, mid=self.maternal_id
        )


def flatten(value, sep=","):
    """
    >>> flatten([1,2,3,4])
    '1,2,3,4'
    >>> flatten((5,6))
    '5,6'
    >>> flatten(0.987654321)
    '0.987654'
    >>> flatten(7)
    '7'
    >>> flatten("flatten")
    'flatten'
    """
    flat = None
    # tuple or list
    if isinstance(value, tuple) or isinstance(value, list):
        flat = sep.join([str(i) for i in value])
    # reformats long float values
    elif isinstance(value, float):
        flat = "%.6f" % (value,)
    # string and int
    else:
        flat = str(value)
    return flat


def zip_lists(value):
    """
    >>> zip_lists([[0,1,2], [3,4,5]])
    ['0 3', '1 4', '2 5']
    """
    return [flatten(i, sep=" ") for i in zip(*value)]


def get_format_fields(ids, variant):
    """
    args:
        ids (list) - list of FORMAT field IDs, e.g. ['AS', 'AP', 'DHFFC']
        variant (pysam.libcbcf.VariantRecord)

    returns:
        list
    """
    fields = list()
    for i in ids:
        fields.append(
            ["%s=%s" % (i, flatten(j.get(i, ""))) for j in variant.samples.values()]
        )
    return zip_lists(fields)


def get_format_title(samples, ids, variant):
    """
    args:
        samples (list) - list of sample IDs in order of VCF annotations
        ids (list) - list of FORMAT field IDs, e.g. ['AS', 'AP', 'DHFFC']
        variant (pysam.libcbcf.VariantRecord)

    returns:
        dict
    """
    fields = get_format_fields(ids, variant)
    return dict(zip(samples, fields))


def make_plot_titles(samples, attr_values):
    """
    keeping this method separate in the event we add more things to the title

    args:
        samples (list) - list of sample IDs
        attr_values (str) - string of VCF FORMAT values

    returns:
        dict

    >>> make_plot_titles(['s1', 's2', 's3'], {'s1': 'AS=0 AP=0', 's2': 'AS=0 AP=1', 's3': 'AS=1 AP=1'})
    {'s1': "'s1 AS=0 AP=0'", 's2': "'s2 AS=0 AP=1'", 's3': "'s3 AS=1 AP=1'"}
    """
    plot_titles = dict()
    for sample in samples:
        if sample in attr_values:
            plot_titles[sample] = quote("%s %s" % (sample, attr_values[sample]))
    return plot_titles


def get_overlap(
    tabix,
    chrom,
    start,
    end,
    priority=["exon", "gene", "transcript", "cds"],
    no_hit="intergenic",
    fix_chr=True,
):
    """
    args:
        tabix (pysam.libctabix.TabixFile) - open TabixFile
        chrom (str)
        start (int)
        end (int)
        priority (Optional[list]) - order of preferred region annotation
        no_hit (Optional[str]) - use this annotation if no matches among priority
        fix_chr (Optional[bool]) - try to fetch a region using both non-'chr' and 'chr' prefix on failures

    returns:
        str
    """
    overlaps = None
    try:
        overlaps = set(
            [i.split("\t")[2].lower() for i in tabix.fetch(chrom, start, end)]
        )
    except IndexError:
        # probably not a gff or gtf
        print("Invalid annotation file specified for --gff")
        overlaps = None
    except ValueError:
        if fix_chr:
            # try removing chr
            if chrom.startswith("chr"):
                overlaps = get_overlap(
                    tabix, chrom[3:], start, end, priority, no_hit, False
                )
            # or adding chr
            else:
                overlaps = get_overlap(
                    tabix,
                    "chr{chrom}".format(chrom=chrom),
                    start,
                    end,
                    priority,
                    no_hit,
                    False,
                )
    except:
        # bad regions
        print(
            "Error fetching {chrom}:{start}-{end}".format(
                chrom=chrom, start=start, end=end
            )
        )
        overlaps = None

    overlap = ""
    if overlaps:
        for feature in priority:
            if feature in overlaps:
                overlap = feature
                break
    else:
        # fetching overlaps failed
        overlap = "unknown"

    if not overlap and no_hit:
        overlap = no_hit
    return overlap


def parse_ped(path, vcf_samples=None):
    if path is None:
        return {}
    samples = []
    look = {}
    for line in open(path):
        samples.append(Sample(line))
        look[samples[-1].id] = samples[-1]

    for s in samples:
        s.dad = look.get(s.paternal_id)
        if s.dad is not None:
            s.dad.kids.append(s)
        s.mom = look.get(s.maternal_id)
        if s.mom is not None:
            s.mom.kids.append(s)
    # match these samples to the ones in the VCF.
    if vcf_samples is not None:
        result = []
        for i, variant_sample in enumerate(vcf_samples):
            if not variant_sample in look:
                continue
            result.append(next(s for s in samples if s.id == variant_sample))
            result[-1].i = i
        samples = result

    return {s.id: s for s in samples}


def get_names_to_bams(bams, name_list=None):
    """
    get mapping from names (read group samples) to bam paths)
    this is useful because the VCF has the names and we'll want the bam paths
    for those samples
    if name_list is passed in as a parameter those will be used instead
    """
    names = {}
    if name_list:
        if len(name_list) != len(bams):
            sys.exit("List of sample IDs does not match list of alignment files.")
        for i, p in enumerate(bams):
            names[name_list[i]] = p
    else:
        for p in bams:
            b = pysam.AlignmentFile(p)
            try:
                names[b.header["RG"][0]["SM"]] = p
            except:
                sys.exit(
                    "No RG field in alignment file "
                    + p
                    + ". \nInclude ordered list of sample IDs to avoid this error"
                )
    return names


def tryfloat(v):
    try:
        return float(v)
    except:
        return v


def to_exprs(astr):
    """
    an expr is just a 3-tuple of "name", fn, value"
    e.g. "DHFFC", operator.lt, 0.7"
    >>> to_exprs("DHFFC < 0.5 & SVTYPE == 'DEL'")
    [('DHFFC', <built-in function lt>, 0.5), ('SVTYPE', <built-in function eq>, 'DEL')]

    >>> to_exprs("CSQ contains 'HIGH'")
    [('CSQ', <built-in function contains>, 'HIGH')]
    """
    astr = (x.strip() for x in astr.strip().split("&"))
    result = []
    for a in astr:
        a = [x.strip() for x in a.split()]
        if len(a) == 2:
            assert a[1] == "exists", ("bad expression", a)
            a.append("extra_arg")
        assert len(a) == 3, ("bad expression", a)
        assert a[1] in cmp_lookup, (
            "comparison:"
            + a[1]
            + " not supported. must be one of:"
            + ",".join(cmp_lookup.keys())
        )
        result.append((a[0], cmp_lookup[a[1]], tryfloat(a[2].strip("'").strip('"'))))
    return result


def check_expr(vdict, expr):
    """
    >>> check_expr({"CSQ": "asdfHIGHasdf"}, to_exprs("CSQ contains 'HIGH'"))
    True

    >>> check_expr({"CSQ": "asdfHIGHasdf", "DHFC": 1.1}, to_exprs("CSQ contains 'HIGH' & DHFC < 0.5"))
    False

    >>> check_expr({"CSQ": "asdfHIGHasdf", "DHFC": 1.1}, to_exprs("CSQ contains 'HIGH' & DHFC < 1.5"))
    True

    >>> check_expr({"smoove_gene": "asdf"}, to_exprs("smoove_gene exists"))
    True

    >>> check_expr({"smooe_gene": "asdf"}, to_exprs("smoove_gene exists"))
    False

    >>> check_expr({"smoove_gene": ""}, to_exprs("smoove_gene exists"))
    True
    """

    # a single set of exprs must be "anded"
    for name, fcmp, val in expr:
        # NOTE: asking for a missing annotation will return false.
        if not name in vdict:
            return False
        if not fcmp(vdict[name], val):
            return False
    return True


def make_single(vdict):
    """
    >>> d = {"xx": (1,)}
    >>> make_single(d)
    {'xx': 1}
    """
    for k in vdict.keys():
        if isinstance(vdict[k], tuple) and len(vdict[k]) == 1:
            vdict[k] = vdict[k][0]
    return vdict


def get_dn_row(ped_samples):
    for s in ped_samples.values():
        if s.mom is not None and s.dad is not None:
            return '{title:"de novo", field:"dn"}'
    return ""


def read_important_regions(bedfilename):
    important_regions = {}
    with open(bedfilename, "r") as bedfile:
        for line in bedfile:
            pos_fields = line.strip().split()
            region_string = "_".join(pos_fields[1:3])
            if pos_fields[0] not in important_regions:
                important_regions[pos_fields[0]] = []
            important_regions[pos_fields[0]].append(region_string)

    return important_regions


def var_in_important_regions(important_regions, chrom, start, end):
    if chrom in important_regions:
        for region in important_regions[chrom]:
            region_st, region_end = [int(x) for x in region.split("_")]
            if (
                region_st <= start <= region_end
                or region_st <= end <= region_end
                or start <= region_st <= end
            ):
                return True
    return False


def cram_input(bams):
    for bam in bams:
        if bam.endswith(".cram"):
            return True
    return False


def vcf(parser, args):
    """
    Generate commands for plotting variants from VCF, create html for showing them
    """

    args, pass_through_args = parser.parse_known_args()
    if args.dn_only and not args.ped:
        sys.exit("Missing --ped, required when using --dn_only")

    if cram_input(args.bams):
        if "-r" not in pass_through_args and not "--reference" in pass_through_args:
            sys.exit(
                "ERROR: missing reference file required for CRAM. "
                + "Use -r option. (Run `samplot.py -h` for more help)"
            )

    global HTML
    global HERE

    vcf = pysam.VariantFile(args.vcf)
    vcf_samples = vcf.header.samples
    vcf_samples_set = set(vcf_samples)
    vcf_samples_list = list(vcf_samples)

    annotations = None
    if args.gff:
        annotations = pysam.TabixFile(args.gff)

    filters = [to_exprs(f) for f in args.filter]

    ped_samples = parse_ped(args.ped, vcf_samples)

    # this is empty unless we have a sample with both parents defined.
    dn_row = get_dn_row(ped_samples)

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    names_to_bams = get_names_to_bams(args.bams, args.sample_ids)
    important_regions = None
    if args.important_regions:
        important_regions = read_important_regions(args.important_regions)
    tabledata = []
    # user requested FORMAT fields to add to plot title
    format_field_ids = None
    if args.format:
        format_field_ids = args.format.split(",")

    out_file = open(args.command_file, "w")

    for variant in vcf:
        svtype = variant.info.get("SVTYPE", "SV")
        if args.important_regions:
            if not var_in_important_regions(
                important_regions, variant.chrom, variant.start, variant.stop
            ):
                continue

        if svtype in ("INS"):
            continue
        if variant.stop - variant.start > args.max_mb * 1000000:
            continue
        if variant.stop - variant.start < args.min_bp:
            continue

        gts = [s.get("GT", (None, None)) for s in variant.samples.values()]

        if sum(None in g for g in gts) >= args.min_call_rate * len(vcf_samples):
            continue
        if args.max_hets:
            # requisite hets/hom-alts
            if sum(sum(x) >= 1 for x in gts if not None in x) > args.max_hets:
                continue
        if not any(sum(x) > 0 for x in gts if not None in x):
            continue

        test_idxs = [i for i, gt in enumerate(gts) if not None in gt and sum(gt) > 0]
        test_samples = [
            s for i, s in enumerate(variant.samples.values()) if i in test_idxs
        ]

        if len(filters) == 0:
            idxs = test_idxs
        else:
            idxs = []
            odict = make_single(dict(variant.info.items()))
            for i, ts in enumerate(test_samples):
                vdict = odict.copy()
                vdict.update(make_single(dict(ts.items())))

                if any(check_expr(vdict, fs) for fs in filters):
                    idxs.append(test_idxs[i])
        if len(idxs) == 0:
            continue
        is_dn = []

        # we call it a de novo if the sample passed the filters but the mom and
        # dad had homref genotypes before filtering.
        # so stringent filtering on the kid and lenient on parents.
        variant_samples = []
        for i in idxs:
            if vcf_samples[i] in names_to_bams:
                variant_samples.append(vcf_samples[i])
        if len(variant_samples) <= 0:
            continue

        bams = [names_to_bams[s] for s in variant_samples]
        if dn_row != "":
            test_sample_names = {s.name for s in test_samples}
            for variant_sample in variant_samples:
                sample = ped_samples[variant_sample]
                if sample.mom is None or sample.dad is None:
                    continue
                if (
                    not sample.mom.id in test_sample_names
                    and not sample.dad.id in test_sample_names
                ):
                    is_dn.append(sample.id)

        if len(is_dn) <= 0 and args.dn_only:
            continue

        # save these for the html.
        n_samples = len(variant_samples)
        # semi-colon delimited eases CSV export from HTML
        sample_str = ";".join(variant_samples)
        # dict holding sample to FORMAT title string
        plot_titles = dict()
        if format_field_ids:
            format_attrs = get_format_title(vcf_samples_list, format_field_ids, variant)
            plot_titles = make_plot_titles(variant_samples, format_attrs)

        # try to get family members
        if args.ped is not None:
            # do DN samples first so we can see parents.
            for variant_sample in is_dn + [
                x for x in variant_samples if not x in is_dn
            ]:
                s = ped_samples.get(variant_sample)
                if s is None:
                    continue
                if (
                    s.mom is not None
                    and not s.mom.id in variant_samples
                    and s.mom.id in vcf_samples_set
                ):
                    variant_samples.append("mom-of-%s[%s]" % (variant_sample, s.mom.id))
                    bams.append(names_to_bams[s.mom.id])
                if (
                    s.dad is not None
                    and not s.dad.id in variant_samples
                    and s.dad.id in vcf_samples_set
                ):
                    variant_samples.append("dad-of-%s[%s]" % (variant_sample, s.dad.id))
                    bams.append(names_to_bams[s.dad.id])
                for kid in s.kids:
                    if not kid.id in variant_samples and kid.id in vcf_samples_set:
                        variant_samples.append(
                            "kid-of-%s[%s]" % (variant_sample, kid.id)
                        )
                        bams.append(names_to_bams[kid.id])
                    if args.max_hets:
                        if len(bams) > 1.5 * args.max_hets:
                            break
                if args.max_hets:
                    if len(bams) > 1.5 * args.max_hets:
                        break
        elif args.min_entries and len(bams) < args.min_entries:
            # extend with some controls:
            hom_ref_idxs = [
                i
                for i, gt in enumerate(gts)
                if len(gt) == 2 and gt[0] == 0 and gt[1] == 0
            ]
            if len(hom_ref_idxs) > 3:
                random.shuffle(hom_ref_idxs)
                hom_ref_idxs = hom_ref_idxs[:3]

            hom_ref_samples = []
            for i in hom_ref_idxs:
                if vcf_samples[i] in names_to_bams:
                    hom_ref_samples.append(vcf_samples[i])

            to_add_count = args.min_entries - len(bams)
            bams.extend(names_to_bams[s] for s in hom_ref_samples[:to_add_count])
            variant_samples += [
                "control-sample:" + s for s in hom_ref_samples[:to_add_count]
            ]

        data_dict = {
            "chrom": variant.chrom,
            "start": variant.start,
            "end": variant.stop,
            "svtype": svtype,
            "svlength": variant.stop - variant.start,
            "samples": sample_str,
            "nsamples": n_samples,
        }
        if annotations:
            data_dict["overlaps"] = get_overlap(
                annotations, variant.chrom, variant.start, variant.stop
            )
        if dn_row != "":
            data_dict["dn"] = ",".join(is_dn)
        fig_path = os.path.join(
            args.out_dir,
            "{svtype}_{chrom}_{start}_{end}.{itype}".format(
                itype=args.output_type, **data_dict
            ),
        )
        tabledata.append(data_dict)

        if "CIPOS" in variant.info:
            v = variant.info["CIPOS"]
            cipos = "--start_ci '%s,%s'" % (abs(v[0]), abs(v[1]))
        else:
            cipos = ""
        if "CIEND" in variant.info:
            v = variant.info["CIEND"]
            ciend = "--end_ci '%s,%s'" % (abs(v[0]), abs(v[1]))
        else:
            ciend = ""
        # dynamically set Z to speed drawing and remove noise for larger events
        z = 3
        if variant.stop - variant.start > 2000:
            z = 4
        if variant.stop - variant.start > 10000:
            z = 6
        if variant.stop - variant.start > 20000:
            z = 9

        if args.max_entries:
            bams = bams[: args.max_entries]
            variant_samples = variant_samples[: args.max_entries]

        # update titles based on FORMAT fields requested
        title_list = list()
        for variant_sample in variant_samples:
            if variant_sample in plot_titles:
                title_list.append(plot_titles[variant_sample])
            else:
                title_list.append(variant_sample)

        out_file.write(
            "samplot plot {extra_args} -z {z} --minq 0 -n {titles} {cipos} {ciend} {svtype} -c {chrom} -s {start} -e {end} -o {fig_path} -d {downsample} -b {bams}\n".format(
                here=HERE,
                extra_args=" ".join(pass_through_args),
                bams=" ".join(bams),
                titles=" ".join(title_list),
                z=z,
                cipos=cipos,
                ciend=ciend,
                svtype="-t " + svtype if svtype != "SV" else "",
                fig_path=fig_path,
                chrom=variant.chrom,
                start=variant.start,
                end=variant.stop,
                downsample=args.downsample,
            )
        )

    if args.command_file:
        out_file.close()

    # update the javascript
    HTML = HTML.replace("[DATA]", json.dumps(tabledata))
    HTML = HTML.replace("[PLOT_TYPE]", args.output_type)
    HTML = HTML.replace("[GFF]", "true" if annotations else "false")
    HTML = HTML.replace("[DENOVO]", "true" if dn_row else "false")

    with open("{out_dir}/index.html".format(out_dir=args.out_dir), "w") as fh:
        print(HTML, file=fh)

    if not args.manual_run:
        import subprocess

        # make runnable and then run it
        os.chmod(args.command_file, 0o755)
        subprocess.call(os.path.join(".", args.command_file), shell=True)
        os.remove(args.command_file)


def add_vcf(parent_parser):
    """Defines allowed arguments for samplot's vcf plotter
    """
    import doctest

    parser = parent_parser.add_parser(
        "vcf",
        help="Generates commands to plot images with `samplot plot`,"
        + " using VCF file to define regions",
    )

    if len(sys.argv) > 1 and sys.argv[1] == "test":
        r = doctest.testmod()
        print(r)
        sys.exit(r.failed)

    parser.add_argument("--vcf", "-v", help="VCF file containing structural variants")
    parser.add_argument(
        "-d", "--out-dir", help="path to write output PNGs", default="samplot-out"
    )
    parser.add_argument("--ped", help="path ped (or .fam) file")
    parser.add_argument(
        "--dn_only",
        help="plots only putative de novo variants (PED file required)",
        action="store_true",
    )
    parser.add_argument(
        "--min_call_rate",
        type=float,
        help="only plot variants with at least this call-rate",
        default=0.95,
    )
    parser.add_argument(
        "--filter",
        action="append",
        help="simple filter that samples"
        + " must meet. Join multiple filters with '&' and specify --filter multiple times for 'or'"
        + " e.g. DHFFC < 0.7 & SVTYPE = 'DEL'",
        default=[],
    )
    parser.add_argument(
        "-O",
        "--output_type",
        choices=("png", "pdf", "eps", "jpg"),
        help="type of output figure",
        default="png",
    )
    parser.add_argument(
        "--max_hets",
        type=int,
        help="only plot variants with at most this many heterozygotes",
        required=False,
    )
    parser.add_argument(
        "--min_entries",
        type=int,
        help="try to include homref samples as controls to get this many samples in plot",
        default=6,
    )
    parser.add_argument(
        "--max_entries",
        type=int,
        help="only plot at most this many heterozygotes",
        default=10,
    )
    parser.add_argument(
        "--max_mb",
        type=int,
        help="skip variants longer than this many megabases",
        default=1,
    )
    parser.add_argument(
        "--min_bp",
        type=int,
        help="skip variants shorter than this many bases",
        default=20,
    )
    parser.add_argument(
        "--important_regions",
        help="only report variants that overlap regions in this bed file",
        required=False,
    )
    parser.add_argument(
        "-b",
        "--bams",
        type=str,
        nargs="+",
        help="Space-delimited list of BAM/CRAM file names",
        required=True,
    )
    parser.add_argument(
        "--sample_ids",
        type=str,
        nargs="+",
        help="Space-delimited list of sample IDs, must have same order as BAM/CRAM file names. BAM RG tag required if this is ommitted.",
        required=False,
    )
    parser.add_argument(
        "--command_file",
        help="store commands in this file.",
        default="samplot_vcf_cmds.tmp",
        required=False,
    )
    parser.add_argument(
        "--format",
        default="AS,AP,DHFFC",
        help="comma separated list of FORMAT fields to include in sample plot title",
    )
    parser.add_argument(
        "--gff",
        help="genomic regions (.gff with .tbi in same directory) used when building HTML table and table filters",
    )
    parser.add_argument(
        "--downsample", help="Number of normal reads/pairs to plot", default=1, type=int
    )
    parser.add_argument(
        "--manual_run",
        help="don't auto-run the samplot plot commands (command_file will be deleted)",
        default=False,
        action="store_true",
    )

    parser.set_defaults(func=vcf)


# if __name__ == "__main__":
#    main(args, pass_through_args)
