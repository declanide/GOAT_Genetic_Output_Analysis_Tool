/**
 * App.js
 */

    // Setup AJAX pour la gestion des token

var csrftoken = $.cookie('csrftoken');

function csrfSafeMethod(method) {
    return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
}
$.ajaxSetup({
    beforeSend: function(xhr, settings) {
        if (!csrfSafeMethod(settings.type) && !this.crossDomain) {
            xhr.setRequestHeader("X-CSRFToken", csrftoken);
        }
    }
});

    //  Ajax sur Submit

//$(document).ready(function() {
//   $("#SGene").submit(function(event){
//        $.ajax({
//             type:"POST",
//             url:"/SGene/",
//             data: {
//                    'rs_ID': $('#rs_ID').val(),
//                    'gene': $('#gene').val(),
//                    'threshold': $('#threshold').val()
//                    },
//             success: function(data){
//                    console.log(data);
//             },
//            error: function(error){
//                console.log(error);
//                $("#error").text(error.responseText);
//            }
//        });
//        return false;
//   });
//});

/*** Form Display ***/

function openForm(name)
{
    id = name;
    closeForm();
    $('#background').fadeIn().css("display","block");
    $('#'+id).addClass('visible');

}

function closeForm()
{
    $('.form').removeClass('visible');
    $('#background').fadeOut().css( "display","none");
}

$( document ).ready(function() {
    $("#background").click(function(){
        closeForm();
    });

    $(".form").click(function(e){
        e.stopPropagation();
    });
});


