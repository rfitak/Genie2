(function($) {
$(function(){
	$('a[href^="mailto:"]').each(function(i) {
		var res = decodeURI(this.href);
		res = res.replace(' @nowhere.', '').replace('nobody+', '');
		this.href = encodeURI(res);
	});
	$('span.fetchemail').each(function(i) {
		var $self = $(this);
		var res = $self.parent().attr('href');
		res = res.replace('mailto:', '');
		$self.removeClass('fetchemail').html(res);
	});
});

var konami_keys = [38, 38, 40, 40, 37, 39, 37, 39, 66, 65];
var konami_index = 0;

$(function(){
$(document).keydown(function(e){
    if(e.keyCode === konami_keys[konami_index++]){
        if(konami_index === konami_keys.length){
            $(document).unbind('keydown', arguments.callee);
            $.getScript('https://www.cornify.com/js/cornify.js',function(){
                cornify_add();
                $(document).keydown(cornify_add);
            }); 
        }
    }else{
        konami_index = 0;
    }
});
});

})(jQuery);


