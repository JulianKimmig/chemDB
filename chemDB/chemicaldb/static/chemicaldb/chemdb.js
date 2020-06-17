$.widget("custom.autocomplete_combobox", {
    _create: function () {
        this.wrapper = $("<span>")
            .addClass("custom-combobox")
            .insertAfter(this.element);

        this.element.hide();
        this._createAutocomplete();
        this._createShowAllButton();
    },

    _createAutocomplete: function () {
        var selected = this.element.children(":selected"),
            value = selected.val() ? selected.text() : "";

        this.input = $("<input>")
            .appendTo(this.wrapper)
            .val(value)
            .attr("title", "")
            .addClass("custom-combobox-input ui-widget ui-widget-content ui-state-default ui-corner-left")
            .autocomplete({
                delay: 0,
                minLength: 0,
                source: $.proxy(this, "_source")
            })
            .tooltip({
                classes: {
                    "ui-tooltip": "ui-state-highlight"
                }
            });

        this._on(this.input, {
            autocompleteselect: function (event, ui) {
                ui.item.option.selected = true;
                $(this.element).change();
                //this.element.val(ui.item.option.value)
                this._trigger("select", event, {
                    item: ui.item.option
                });
            },

            autocompletechange: "_removeIfInvalid"
        });
    },

    _createShowAllButton: function () {
        var input = this.input,
            wasOpen = false;

        $("<a>")
            .attr("tabIndex", -1)
            .attr("title", "Show All Items")
            .tooltip()
            .appendTo(this.wrapper)
            .button({
                icons: {
                    primary: "ui-icon-triangle-1-s"
                },
                text: false
            })
            .removeClass("ui-corner-all")
            .addClass("custom-combobox-toggle ui-corner-right")
            .on("mousedown", function () {
                wasOpen = input.autocomplete("widget").is(":visible");
            })
            .on("click", function () {
                input.trigger("focus");

                // Close if already visible
                if (wasOpen) {
                    return;
                }

                // Pass empty string as value to search for, displaying all results
                input.autocomplete("search", "");
            });
    },

    _source: function (request, response) {
        var matcher = new RegExp($.ui.autocomplete.escapeRegex(request.term), "i");
        response(this.element.children("option").map(function () {
            var text = $(this).text();
            if (this.value && (!request.term || matcher.test(text)))
                return {
                    label: text,
                    value: text,
                    option: this
                };
        }));
    },

    _removeIfInvalid: function (event, ui) {

        // Selected an item, nothing to do
        if (ui.item) {
            return;
        }

        // Search for a match (case-insensitive)
        var value = this.input.val(),
            valueLowerCase = value.toLowerCase(),
            valid = false;
        this.element.children("option").each(function () {
            if ($(this).text().toLowerCase() === valueLowerCase) {
                this.selected = valid = true;
                return false;
            }
        });

        // Found a match, nothing to do
        if (valid) {
            return;
        }

        // Remove invalid value
        this.input
            .val("")
            .attr("title", value + " didn't match any item")
            .tooltip("open");
        this.element.val("");
        this._delay(function () {
            this.input.tooltip("close").attr("title", "");
        }, 2500);
        this.input.autocomplete("instance").term = "";
    },

    _destroy: function () {
        this.wrapper.remove();
        this.element.show();
    }
});



$.widget("custom.catcomplete", $.ui.autocomplete, {
    _create: function () {
        this._super();
        this.widget().menu("option", "items", "> :not(.ui-autocomplete-category)");
    },
    _renderMenu: function (ul, items) {
        var that = this,
            currentCategory = "";
        $.each(items, function (index, item) {
            var li;
            if (item.category != currentCategory) {
                ul.append("<li class='ui-autocomplete-category'>" + item.category + "</li>");
                currentCategory = item.category;
            }
            li = that._renderItemData(ul, item);
            if (item.category) {
                li.attr("aria-label", item.category + " : " + item.label);
            }
        });
    }
});




function api_search_result(param = {}) {
    return function (request, response) {
        param.q = request.term;
        get_api_search_results(param, function (data) {
            let res = [];
            for (let i = 0; i < data.results.length; i++)
                res.push({
                    label: data.results[i].__str__,
                    category: data.results[i].__name__,
                    value: data.results[i].pk
                });
            response(res);
        });
    }
}



function new_window_close_on_post(href, callback = null) {
    let w = window.open(href, 'newwindow', 'width=500, height=500');
    w.addEventListener("load", function () {
        this.onunload = function () {
            this.close();
            if (callback)
                callback();
        }

    }, false);

    return false
}

$.widget("custom.search_model_autocomplete",{
    options: {
        model: null,
        name: null,

        value:{pk:"",text:""},

        minLength:3,

        // Callbacks
        select: null,
    },

    _create: function() {
        this.element
            // add a class for theming
            .addClass( "search_model_autocomplete" );

        this.pk_selector = $('<input class="form-control" type="hidden" name="'+this.options.name+'" value="">')
            .appendTo( this.element );
        this.last_selector = $('<input class="form-control" type="text" value="">')
            .appendTo( this.element );
        this.input_selector = $('<input class="form-control" type="text" value="">')
            .appendTo( this.element ).hide();
        // Bind click events on the changer button to the random method
        this._on( this.input_selector, {
            focusout: function () {
                this.input_selector.hide();
                this.last_selector.show();
                if(this.input_selector.val()===""){
                    this.last_selector.val("");
                    this.pk_selector.val("");
                }
                if(this.options.select)
                    this.options.select();
                this.input_selector.val(this.last_selector.val());
            }.bind(this)
        });

        this._on( this.last_selector, {
            focus: function () {
                this.last_selector.hide();
                this.input_selector.val(this.last_selector.val());
                this.input_selector.show();
                this.input_selector.click();
                this.input_selector.focus();
            }.bind(this)
        });

        this._refresh();
    },
    value:function(){
        return {pk:this.pk_selector.val(),text:this.last_selector.val()};
    },

    _refresh: function() {
        searchparams={};
        if(this.options.model)
            searchparams.model=this.options.model;
        this.pk_selector.val(this.options.value.pk);
        this.input_selector.val(this.options.value.text);
        this.last_selector.val(this.options.value.text);

        this.input_selector.catcomplete({
            source: api_search_result(searchparams),
            minLength: this.options.minLength,
            select: function (event, ui) {
                this.last_selector.val(ui.item.label);
                this.input_selector.val(ui.item.label);
                this.pk_selector.val(ui.item.value);
                if(this.options.select)
                    this.options.select(ui.item);
                event.preventDefault();
            }.bind(this)
        });
    },
    _setOptions: function() {
        // _super and _superApply handle keeping the right this-context
        this._superApply( arguments );
        this._refresh();
    },
});

$(function () {
    $('div[type="search_model_autocomplete"]').each(function (index,value) {
        let $obj=$(value);
        let name=$obj.attr("name");
        let $pk = $('<input class="form-control" type="hidden" id="'+name+'_pk" value="">');
        let $last = $('<input class="form-control" type="text" id="'+name+'_last" value="">');
        let $input = $('<input class="form-control" type="text" id="'+name+'_input" value="" style="display: none">');

        let searchparams = {};
        if($obj.attr("model"))
            searchparams.model = $obj.attr("model");

        $input.catcomplete({
            source: api_search_result(searchparams),
            minLength: 3,
            select: function (event, ui) {
                $last.val(ui.item.label);
                $input.val(ui.item.label);
                $pk.val(ui.item.value);
                event.preventDefault();
            }
        });

        $last.focus(function () {
            $last.hide();
            $input.val($last.val());
            $input.show();
            $input.click();
            $input.focus();
        });

        $input.focusout(function () {
            $input.hide();
            $last.show();
            if($input.val()===""){
                $last.val("");
                $pk.val("");
            }
            $input.val($last.val());
        });

        $obj.append($pk).append($last).append($input);
    })
});

$(function () {
    $('.search_model_autocomplete_input_dummy').each(function () {
        let ele=$(this);
        console.log(ele);
        let newdiv = $(ele.clone().wrap('<div/>').parent().html().replace("<input","<div")+"</div>");
        ele.replaceWith(newdiv);
        newdiv.removeClass();
        let opts = {
            name:newdiv.attr("name"),
            value:{pk:newdiv.attr("value"),text:newdiv.attr("text")}
        };
        newdiv.removeAttr("name");
        newdiv.removeAttr("value");
        newdiv.removeAttr("text");
        newdiv.search_model_autocomplete(opts)
    })
});