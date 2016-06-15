page = '';

$(function() {
    $('#password').modal({
        backdrop: 'static',
        keyboard: false
    });
    var pathname = window.location.pathname;
    page = pathname.split('/')[2];
    $.get('/' + page, function(data, status) {
        if (status == 'success') {
            $("#editor").val(data);
        }
    });
    $('#enter').click(function() {
        $.get('/check_pwd/' + $('#pwd').val(), function(data, status) {
            if (data == 'success') {
                $('#password').modal('hide');
            } else {
                alert("wrong password");
            }
        });

    });

    $("#save").click(function() {
        $.post('/save/' + page, { data: $('#editor').val() }, function(data, status) {
            if (data == 'success') {
                alert("page has been saved!");
            } else {
                alert("failed to save page!");
            }
        });
    });

});
