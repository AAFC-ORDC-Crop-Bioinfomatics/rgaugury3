prj_id = ""
gene_img_1 = '<a href="'
gene_img_2 = 'gallery/{data}" target="_blank"><img src ="'
gene_img_3 = '{data}" width="100%" height="100%"></a>'
invalid_img = '<a href="#" class="thumbnail disabled"> <span class="glyphicon glyphicon-picture"></span> </a>'
table_header = undefined
nbs_data = undefined
rlp_data = undefined
rlk_data = undefined
tmcc_data = undefined
root=""

$(function() {
    prj_id = $("#prj_id").text().trim();
    root = $("#root").text().trim();

    highlight($('#gene_type button').first());
    
    $('#gene_type button').click(function() {
        highlight($(this))
        initDSTable(activeTable, $(this).text());
      
    });

//    $(document).on('click', '.ds_table img', function(event) {
//        var src = $(this).attr('src');
//        var url = '/gallery' + src;
//        window.location.href = url;
//        event.preventDefault() 
//    });

    $('#gallery').bind('DOMMouseScroll mousewheel', function(event) {
        var gallery = $('#gallery .modal-dialog');
        var image = $('.modal-body img');
        var width = gallery.width();
        if (event.originalEvent.wheelDelta > 0 || event.originalEvent.detail < 0) {
            // scroll up
            gallery.width(width * 1.1);
        } else {
            // scroll down
            gallery.width(width * 0.9);
        }
    });

    initDSTable(getAlldata);
    $('#domain th:nth-child(4)').click(function(event) {
        event.stopPropagation();
    });
});

function highlight(button) {
    var siblings = $('#gene_type button').not(button);
    button.css('background-color', '#337ab7');
    button.css('color', 'white');
    siblings.css('background-color', 'white');
    siblings.css('color', '#337ab7');
}

function getAlldata() {
    getTableData('NBS', true);
    getTableData('RLP');
    getTableData('RLK');
    getTableData('TM-CC');
}

function initDSTable(callback, type) {
    var domain = $('#domain');
    domain.html('');
    domain.html('<table id="ds_table" class="table table-striped table-bordered ds_table" width="100%" cellspacing="0"></table>');
    if (table_header) {
        $('#ds_table').append(table_header);
        callback(type);
    } else {
        addHeader(callback);
    }
}

function activeTable(type) {
    var dataType = type.toLowerCase().replace(/-/, '');
    var data = undefined;
    switch (dataType.trim()) {
        default:
            case 'nbs':
            data = nbs_data;
        break;
        case 'rlp':
                data = rlp_data;
            break;
        case 'rlk':
                data = rlk_data;
            break;
        case 'tmcc':
                data = tmcc_data;
            break;
    }

    setDataTable(data);
}

function addHeader(callback) {
    $.ajax({
        "url": root+'ds_header/' + prj_id,
        "dataType": "json",
        "success": function(json) {
            var tableHeaders = '';
            $.each(json.columns, function(i, val) {
                tableHeaders += '<th width="150px">' + val.trim() + "</th>";
            });
            table_header = '<thead><tr>' + tableHeaders + '</tr></thead>';
            $('#ds_table').append(table_header);
            $('#ds_table th:eq(3)').click(function(e) {
                e.stopImmediatePropagation();
            });
            callback();
        }
    });
}
// if setTable is true, then initialize the table with the json value
function getTableData(type, setTable) {
    $.ajax({
        "url": root+'gene_info/' + prj_id + '/' + type,
        "dataType": "json",
        "success": function(json) {
            var dataType = type.toLowerCase().replace(/-/, '');
            switch (dataType) {
                default:
                    case 'nbs':
                    nbs_data = json.data;
                break;
                case 'rlp':
                        rlp_data = json.data;
                    break;
                case 'rlk':
                        rlk_data = json.data;
                    break;
                case 'tmcc':
                        tmcc_data = json.data;
                    break;
            }
            if (setTable) {
                setDataTable(json.data);
            }
        }
    });
};

function addXSclass() {
    $('.col-sm-5').addClass('col-xs-5')
    $('.col-sm-6').addClass('col-xs-6')
    $('.col-sm-7').addClass('col-xs-7')
}

function setDataTable(data) {
    $('#ds_table').DataTable({
        "data": data,
        "createdRow": function(row, data, index) {
            var colum = 3;
            var img_path = data[colum];
            if (img_path.search(/img/i) < 0) {
                $('td', row).eq(colum).html(invalid_img);
            } else {
                var src = 'img/' + prj_id + '/' + data[colum];
                $('td', row).eq(colum).html(gene_img_1+root+gene_img_2.replace(/{data}/g, src)+root+gene_img_3.replace(/{data}/g, src));
            }
        },
        "initComplete": function(settings, json) {
            locateFooter();
            addXSclass();
        },
        "fnDrawCallback": function(oSettings) {
            locateFooter();
        },
        "pagingType": "input",
        "aaSorting": [[2,'asc'],[0,'asc']],
    });
}
