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

$(function () {
    $('div[type="search_model_autocomplete"]').each(function (index,value) {
        let $obj=$(value);
        let name=$obj.attr("name");
        let $pk = $('<input class="form-control" type="hidden" id="'+name+'_pk" value="">');
        let $last = $('<input class="form-control" type="text" id="'+name+'_last" value="">');
        let $input = $('<input class="form-control" type="text" id="'+name+'_input" value="" style="display: none">');

        $input.catcomplete({
            source: api_search_result({model: "Substance"}),
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

        });

        $input.focusout(function () {
            $input.hide();
            $last.show();
            $input.val($last.val());
        });

        $obj.append($pk).append($last).append($input);
    })
});