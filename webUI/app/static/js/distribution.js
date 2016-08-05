$(function(){
	$('#list li').click(function(){
		var text =$(this).find('div').first().text();
		var src = $(this).find('img').attr('src');
		$('#type h3').text(text);
		$('#img img').attr('src',src);
	});
    $('#img img').dblclick(function(){
         window.open($(this).attr('src'), '_blank');
    });
});