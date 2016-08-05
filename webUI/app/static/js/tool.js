function download_to_textbox(tringer, url, textarea, sucess_func, fail_func) {
    $.get(url, null, function(data, status) {
        if (status == 'success') {
            textarea.val(data);
            sucess_func();
        } else {
            fail_func();
        }
    }, "text");
}

String.prototype.hashCode = function(){
 var hash = 0;
 if (this.length == 0) return hash;
 for (i = 0; i < this.length; i++) {
  char = this.charCodeAt(i);
  hash = ((hash<<5)-hash)+char;
  hash = hash & hash; // Convert to 32bit integer
 }
 return hash;
}

function fingerprint() {
    var canvas = document.createElement('canvas');
    var ctx = canvas.getContext('2d');
    var txt = 'i9asdm..$#po((^@KbXrww!~cz';
    ctx.textBaseline = "top";
    ctx.font = "16px 'Arial'";
    ctx.textBaseline = "alphabetic";
    ctx.rotate(.05);
    ctx.fillStyle = "#f60";
    ctx.fillRect(125,1,62,20);
    ctx.fillStyle = "#069";
    ctx.fillText(txt, 2, 15);
    ctx.fillStyle = "rgba(102, 200, 0, 0.7)";
    ctx.fillText(txt, 4, 17);
    ctx.shadowBlur=10;
    ctx.shadowColor="blue";
    ctx.fillRect(-20,10,234,5);
    var strng=canvas.toDataURL();
    return strng.hashCode();
}
