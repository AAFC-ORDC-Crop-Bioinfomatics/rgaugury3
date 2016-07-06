vars = {};
OP_DEL = 'delete';
OP_CANCEL = 'cancel';
root=""

function onStart() {
    root = $("#root").text().trim();

    vars.fingerprint = fingerprint();

    $(document).on('click', '.cancel', function() {
        if ($(this).hasClass('disabled')) {
            return;
        }
        vars.prj_id = $(this).closest('tr').find('td:eq(0)').text().trim();
        vars.btn = $(this);
        $.get(root+'fingerprint/' + vars.prj_id + '/' + vars.fingerprint, function(data, status) {
            if (status == 'success') {
                if (data == '0') {
                    // finger print match
                    var prj_name = $(this).closest('tr').find('td:eq(1)').text().trim();
                    vars.operation = OP_CANCEL;
                    $('span#prj_name').text(prj_name);
                    $('span#op_name').text(OP_CANCEL);
                    $('#modal').modal();
                } else {
                    // finger print does not match
                    $('#fingerprint_modal').modal('show');
                }
            } else {
                alert("failed to check fingerprint");
            }
        });
    });

    $(document).on('click', '.del_prj', function() {
        if ($(this).hasClass('disabled')) {
            return;
        }
        vars.prj_id = $(this).closest('tr').find('td:eq(0)').text().trim();
        vars.btn = $(this);
        $.get(root+'fingerprint/' + vars.prj_id + '/' + vars.fingerprint, function(data, status) {
            if (status == 'success') {
                if (data == '0') {
                    // finger print match
                    var prj_name = $(this).closest('tr').find('td:eq(1)').text();
                    vars.operation = OP_DEL;
                    $('span#prj_name').text(prj_name);
                    $('span#op_name').text(OP_DEL);
                    $('#modal').modal();
                } else {
                    // finger print does not match
                    $('#fingerprint_modal').modal('show');
                }
            } else {
                alert("failed to check fingerprint");
            }
        });
    });

    $(document).on('click', 'button#confirm', function() {
        $('#modal').modal('hide');
        vars.btn.addClass('disabled');
        vars.btn.text('wait...');
        if (vars.operation == OP_CANCEL) {
            $.post("cancel", {
                prj_id: vars.prj_id
            }, function(data, status) {
                if (status == 'success') {
                    vars.btn.parent().closest('tr').find('button.del_prj').removeClass('disabled');
                    td_status = vars.btn.closest('tr').find('td:eq(2)');
                    td_status.text('canceled');
                    td_status.css('color', 'red');
                    vars.btn.text('Cancel');
                }
            });
        } else if (vars.operation == OP_DEL) {
            $.post("delete", {
                prj_id: vars.prj_id
            }, function(data, status) {
                if (status == 'success') {
                    var delay = 1000;
                    vars.btn.closest('tr').fadeOut(delay).promise().done(locateFooter);
                }
            });
        }
    });

    $(document).on({
        mouseenter: function() {
            //stuff to do on mouse enter
        },
        mouseleave: function() {
            //stuff to do on mouse leave
            $('.tooltip').remove();
        }
    }, ".progress");

    $('body').tooltip({
        'selector': '.progress',
        'placement': 'top',
        'html': true,
        'container': 'body'
    });

//    setInterval(update, 1000);
}

$(onStart);

function update() {
    $.get('/all-prj', function(data, status) {
        if (status == 'success') {
            $('#prj_table').html(data);
            if (!$('.progress').length) {
                $('.tooltip').remove();
            }
        }
    });
}
