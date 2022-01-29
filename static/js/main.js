(function($) {

    var form = $("#signup-form");
    jQuery.validator.addMethod("nucleotide_only", function(value, element) {
      return this.optional(element) || /^[acgtACGT]*$/i.test(value);
      }, "Nucleotide only please");
    
      
      form.validate (  {
        // Specify validation rules
        rules: {
          // The key name on the left side is the name attribute
          // of an input field. Validation rules are defined
          // on the right side
          
      
      
          Sequence: { required:true, minlength:20, maxlength:30, nucleotide_only:true    } ,
          
        },
        // Specify validation error messages
        messages: {
          Sequence: { nucleotide_only:"Nucleotide only please", minlength:"Enter sequence 20bp upto 30bp of length",maxlength:"Enter sequence 20bp upto 30bp of length" },
            // Make sure the form is submitted to the destination defined
        // in the "action" attribute of the form when valid
      
      
        onFinished: function(event, currentIndex) {
                alert('Submited');  
                                                  }
                                                }


        });          




    form.steps({
        headerTag: "h3",
        bodyTag: "fieldset",
        transitionEffect: "fade",
        autoFocus: "true",

        labels: {
            previous: 'Previous',
            next: 'Next',
            finish: 'Finish',
            current: ''
        },
        titleTemplate: '<h3 class="title">#title#</h3>'
          });





    
//    $(".toggle-password").on('click', function() {

  //      $(this).toggleClass("zmdi-eye zmdi-eye-off");
    //    var input = $($(this).attr("toggle"));
      //  if (input.attr("type") == "password") {
        //  input.attr("type", "text");
        //} else {
          //input.attr("type", "password");
        //}
    //});

})(jQuery);