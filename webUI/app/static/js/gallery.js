section_info = undefined
prj_id = ''
gene_name = ''
scale = 0
data = undefined
root=""

img_tag = '<img src ="{src}" alt="{src}" />'

$(window).on("load", function() {
    root = $("#root").text().trim();
    prj_id = $("#prj_id").text().trim();
    gene_name = $("#gene_name").text().trim();

    initTable();
    $('a').click(function(event){
        location.reload()
    });

    gene = $('.gene_img');
    var height = gene.prop('naturalHeight');
    var width = gene.prop('naturalWidth');

    drawRuler(gene.width(), width);

    scale = 1.0 * width / gene.width();
    var scaled = '1';

    $('#scale').text(scaled + ':' + scale);

    $.get(root+'section/' + prj_id + '/' + gene_name, function(data, status) {
        if (status == 'success') {
            window.data = data.data
        }
    });

    gene.mouseleave(function() {
        $('.hinttext').css('visibility', 'hidden');
    });

    $(window).resize(function() {

    });

    gene.mousemove(function(event) {
        position = event.pageX - $(this).offset().left;
        n_position = scale * position;

        if (data == undefined) {
            var hint = $('.hinttext');
            $('#type').text('Loading...');
            hint.css('top', -hint.height());
            hint.css('left', position);
            hint.css('visibility', 'visible');
            return;
        }

        info = getInfo(n_position);
        var hint = $('.hinttext');
        if (info) {
            $('#type').text(info[2]);
            $('#span').text(info[0] + ' - ' + info[1]);
            hint.css('top', -hint.height());
            hint.css('left', position);
            hint.css('visibility', 'visible');
        } else {
            hint.css('visibility', 'hidden');
        }
    });
});


function getInfo(position) {
    if (data == undefined) {
        return undefined;
    }
    for (var i in data) {
        var list = data[i][0].split('|');
        var min = list[0] | 0;
        var max = list[1] | 0;
        if (position >= min && position <= max) {
            return [min, max, data[i][1]];
        }
    }
}

function initTable() {
    $('#section_table').DataTable({
        "aaSorting": [0, 'asc'],
        "createdRow": function(row, data, index) {
            var colum = 2;
            var img_path = data[colum];
            $('td', row).eq(colum).html(img_tag.replace(/{src}/g, img_path));
        },
        "pageLength": 25
    });

    $('#gff_table').DataTable({
        "aaSorting": [0, 'asc'],
        "pageLength": 10
    });
}

//width is the whole actual pixel
//max is the max value dispaled on the ruler
function drawRuler(width, max) {
    // configure
    var height = 8
    var font_size = 15
    var max_block = 100
    var each_1_height = 3
    var each_5_height = 5
    var each_10_height = 7
    var offset = 0
        // end configure

    // init
    var org_max = max;
    var org_width = width;
    var remain = 0;
    var new_max = 0;
    if ((remain = max % 1000) != 0) {
        new_max = max - remain + 1000;
    }
    width = new_max * width / max;
    max = new_max
    var canvas = document.getElementById('ruler');
    var context = canvas.getContext('2d');
    context.font = font_size + 'px Calibri';
    canvas.width = org_width + font_size * 4;
    canvas.height = 30;


    // draw main line
    //context.beginPath();
    //context.moveTo(0, height);
    //context.lineTo(width, height);
    //context.stroke();

    // draw blocks
    var span = width / max_block;
    var span_value = max / max_block;
    var lastWidth = 0;
    for (var i = 0; i <= max_block; i++) {
        var distance = span * i;
        var distance_value = span_value * i;
        lastWidth = distance;
        if (org_max < distance_value) {
            break;
        }
        context.beginPath();

        if (i % 10 == 0) {
            block_hight = each_10_height;
            context.fillText(distance_value, distance - offset, font_size + height);
        } else if (i % 5 == 0) {
            block_hight = each_5_height;
            context.fillText(distance_value, distance, font_size + height);
        } else {
            block_hight = each_1_height;
        }
        var calcedHeight = height + block_hight / 2;
        context.moveTo(distance, calcedHeight);
        context.lineTo(distance, calcedHeight - block_hight);
        context.stroke();
    }
    context.fillText('(bp)', lastWidth + font_size, font_size + height);

}
