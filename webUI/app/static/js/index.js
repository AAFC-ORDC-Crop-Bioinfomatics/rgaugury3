ips_url = 'https://www.ebi.ac.uk/interpro/interproscan.html';
versionIsOk = false;
versionStatus = 'Please wait, Checking the latest version of Interproscan...';
STILL_RUN = 'Your interproscan is not the latest verion or there is no the Internet connection. Would you like to continue any way? ';
v_status = 0;
count = 0;
$(function() {
    init();
    calcFingerPrint();
    checkversion();

    function success() {
        $('.sample_fasta').text('Reload sample');
        $('#seq_file').val('');
        $('.sample_fasta').removeClass('disabled');
    }

    function fail() {
        $('.sample_fasta').text('Error!');
    }

    $('.sample_fasta').click(function() {
        $(this).text('Loading ...');
        $(this).addClass('disabled');
        $('#protein_seq').removeClass('place-holder');
        download_to_textbox($(this), "/sample_fasta", $('#protein_seq'), success, fail);
        sampleGff(true)
    });

    $('#seq_file').change(function(e) {
        if (e.target.value === '') {
            // do nothing
        } else {
            $("#protein_seq").val('');
            sampleGff(false)
        }
    });

    $('#gff3_file').change(function(e) {
        if (e.target.value === '') {
            // do nothing
        } else {
            sampleGff(false)
        }
    });


    $(window).click(function(e) {
        if ($('#protein_seq').val() == '') {
            if (!$('#protein_seq').is(":focus")) {

            }
        }
    });


    $('#protein_seq').on('input', function(e) {
        if (e.target.value === '') {
            // Textarea has no value
            // do nothing
        } else {
            // Textarea has a value
            $("#seq_file").val('');
            sampleGff(false)
        }
    });

    $('#protein_seq').focus(function(e) {
        if ($(this).hasClass('place-holder')) {
            $(this).removeClass('place-holder');
            emtpyBlack($(this));
        }
    });
    $('#protein_seq').focusout(function(e) {
        if ($(this).val() == '') {
            loadSampleFasta();
        }
    });

    $('#ev').focus(function(e) {
        if ($(this).hasClass('place-holder')) {
            $(this).removeClass('place-holder');
            emtpyBlack($(this));
        }
    });
    $('#ev').focusout(function(e) {
        if ($(this).val() == '') {
            loadDefaultEvalue();
        }
    });

    $('#proj_name').focus(function(e) {
        if ($(this).hasClass('place-holder')) {
            $(this).removeClass('place-holder');
            emtpyBlack($(this));
        }
    });
    $('#proj_name').focusout(function(e) {
        if ($(this).val() == '') {
            loadDefaultPrjName();
        }
    });

    $('#still-run').click(function() {
        versionIsOk = true;
        $('form').trigger('submit');
    });

    $('button[type="submit"]').click(function(e) {
        if (versionIsOk) {
            // the version of interproscan is the latest
        } else {
            e.preventDefault();
            if (v_status == -1) {
                $('#modal h4').text(STILL_RUN);
                $('#still-run').css('visibility', 'visible');
                $('#modal').modal('show');
            } else {
                $('#modal h4').text(versionStatus);
                $('#modal').modal('show');
                checkversion();
            }

        }
    });

    $('button[type="reset"]').click(function(e) {
        loadDefaultPrjName();
        loadDefaultEvalue();
        loadSampleFasta();
        e.preventDefault();
    });

    handleResize($('textarea'));

    $('input.disabled').parent().css("color", "#aaaaaa");
});

function handleResize(textarea) {
    var resizeInt = null;

    // the handler function
    var resizeEvent = function() {
        locateFooter();
    };

    textarea.on('mousedown', function(e) {
        resizeInt = setInterval(resizeEvent, 100);

    });
    $(window).on("mouseup", function(e) {
        if (resizeInt !== null) {
            clearInterval(resizeInt);
        }
        resizeEvent();
    });
}
// get latest interproscan version
function checkversion() {
    $.ajax({
        url: '/latestVersion',
        dataType: "text",
        success: function(data) {
            if (data == '1') {
                // need to upgrade interproscan
                versionIsOk = false;
                versionStatus = 'Please upgrade interproscan first.'
                v_status = -1;
            } else if (data == '-1') {
                versionIsOk = false;
                versionStatus = 'Failed to check latest interproscan version';
                v_status = -1;
            } else if (data == '0') {
                versionIsOk = true;
                if ($('#modal').hasClass('in')){
                   $('form').trigger('submit'); 
                } 
            }
        }
    });
}

function calcFingerPrint() {
    $('input[name="fingerprint"]').val(fingerprint());
}

function init() {
    loadSampleFasta();
    loadDefaultEvalue();
    loadDefaultPrjName();
}

function loadDefaultPrjName() {
    $('#proj_name').addClass('place-holder');
    $('#proj_name').val('RGA');
}

function loadDefaultEvalue() {
    $('#ev').addClass('place-holder');
    $('#ev').val('1e-5');
}

function loadSampleFasta() {
    $.get('/sample_fasta', function(data, status) {
        if (status = 'success') {
            $('#protein_seq').addClass('place-holder');
            $('#protein_seq').val(data);
        }
    });
    sampleGff(true);
}

function emtpyBlack(element) {
    element.val('');
    element.css('color', 'black');
}

function sampleGff(enabled){
    if (enabled){
        $('#div-sample-gff').show()
        $('#div-gff').hide()
        $("input[name='gff3']").val('default')
    }else{
        $('#div-sample-gff').hide()
        $('#div-gff').show()
        $("input[name='gff3']").val("")
    }
}