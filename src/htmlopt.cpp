#include "htmlopt.h"

void HtmlOpt::init(const std::string& cfg){
    if(cfg.empty() || (!util::exists(cfg))){// initialize use internet cdn
        // js
        hashjs[schjs.bootstrap_bundle_js] = SRC_BOOTSTRAP_BUNDLE_JS;
        hashjs[schjs.jszip_js] = SRC_JSZIP_JS;
        hashjs[schjs.plotly_js] =  SRC_PLOTLY_JS;
        hashjs[schjs.jquery_js] = SRC_JQUERY_JS;
        hashjs[schjs.buttons_html5_js] = SRC_BUTTONS_HTML5_JS;
        hashjs[schjs.jquery_dataTables_js] = SRC_JQUERY_DATATABLES_JS;
        hashjs[schjs.dataTables_buttons_js] = SRC_DATABLES_BUTTONS_JS;
        hashjs[schjs.logo_bundle_js] = SRC_LOGO_JS;
        hashjs[schjs.druid_js] = SRC_DRUID_JS;
        // css
        hashcss[schcss.bootstrp_css] = SRC_BOOTSTRAP_MIN_CSS;
        hashcss[schcss.jquery_dataTables_css] = SRC_QUERY_DATATABLES_CSS;
        hashcss[schcss.buttons_dataTables_css] = SRC_DATABLES_BUTTONS_CSS;
    }else{// use configured cdn
        // read json
        int len = 0, max = 0, tmp;
        char* json = NULL, buf[0x10000];
        FILE* fp = fopen(cfg.c_str(), "rb");
        while((tmp = fread(buf, 1, 0x10000, fp)) != 0){
            if(len + tmp + 1 > max){
                  max = len + tmp + 1;
                  kroundup32(max);
                  json = (char*)realloc(json, max);
              }
              memcpy(json+len, buf, tmp);
              len += tmp;
        }
        fclose(fp);
        kson_t* kson = kson_parse(json);
        free(json);
        // begin parse
        const kson_node_s* rn = kson_by_path(kson->root, 1, schjs.bootstrap_bundle_js.c_str());
        hashjs[schjs.bootstrap_bundle_js] = rn->v.str;
        rn = kson_by_path(kson->root, 1, schjs.jszip_js.c_str());
        hashjs[schjs.jszip_js] = rn->v.str;
        rn = kson_by_path(kson->root, 1, schjs.plotly_js.c_str());
        hashjs[schjs.plotly_js] = rn->v.str;
        rn = kson_by_path(kson->root, 1, schjs.jquery_js.c_str());
        hashjs[schjs.jquery_js] = rn->v.str;
        rn = kson_by_path(kson->root, 1, schjs.buttons_html5_js.c_str());
        hashjs[schjs.buttons_html5_js] = rn->v.str;
        rn = kson_by_path(kson->root, 1, schjs.jquery_dataTables_js.c_str());
        hashjs[schjs.jquery_dataTables_js] = rn->v.str;
        rn = kson_by_path(kson->root, 1, schjs.dataTables_buttons_js.c_str());
        hashjs[schjs.dataTables_buttons_js] = rn->v.str;
        rn = kson_by_path(kson->root, 1, schjs.logo_bundle_js.c_str());
        hashjs[schjs.logo_bundle_js] = rn->v.str;
        rn = kson_by_path(kson->root, 1, schjs.druid_js.c_str());
        hashjs[schjs.druid_js] = rn->v.str;
        rn = kson_by_path(kson->root, 1, schcss.bootstrp_css.c_str());
        hashcss[schcss.bootstrp_css] = rn->v.str;
        rn = kson_by_path(kson->root, 1, schcss.jquery_dataTables_css.c_str());
        hashcss[schcss.jquery_dataTables_css] = rn->v.str;
        rn = kson_by_path(kson->root, 1, schcss.buttons_dataTables_css.c_str());
        hashcss[schcss.buttons_dataTables_css] = rn->v.str;
        // destroy
        kson_destroy(kson);
    }
}

void HtmlOpt::printJsAndCSSRes(kstring_t* s){
    // css
    ksprintf(s, "<link type=\"text/css\" href=\"%s\" rel=\"stylesheet\">\n", hashcss[schcss.bootstrp_css].c_str());
    ksprintf(s, "<link type=\"text/css\" href=\"%s\" rel=\"stylesheet\">\n", hashcss[schcss.buttons_dataTables_css].c_str());
    ksprintf(s, "<link type=\"text/css\" href=\"%s\" rel=\"stylesheet\">\n", hashcss[schcss.jquery_dataTables_css].c_str());
    // js
    ksprintf(s, "<script src='%s'></script>\n", hashjs[schjs.jquery_js].c_str());
    ksprintf(s, "<script src='%s'></script>\n", hashjs[schjs.bootstrap_bundle_js].c_str());
    ksprintf(s, "<script src='%s'></script>\n", hashjs[schjs.jquery_dataTables_js].c_str());
    ksprintf(s, "<script src='%s'></script>\n", hashjs[schjs.dataTables_buttons_js].c_str());
    ksprintf(s, "<script src='%s'></script>\n", hashjs[schjs.jszip_js].c_str());
    ksprintf(s, "<script src='%s'></script>\n", hashjs[schjs.buttons_html5_js].c_str());
    ksprintf(s, "<script src='%s'></script>\n", hashjs[schjs.plotly_js].c_str());
    ksprintf(s, "<script src='%s'></script>\n", hashjs[schjs.logo_bundle_js].c_str());
    ksprintf(s, "<script src='%s'></script>\n", hashjs[schjs.druid_js].c_str());
}

void HtmlOpt::printExportButtons(kstring_t* s, std::string tid, std::string fname){
    ksprintf(s, "<script>\n");
    ksprintf(s, "    var rowCount = $('#%s >tbody >tr').length;\n", tid.c_str());
    ksprintf(s, "    var lengthMenuAll = [[rowCount], [\"All\"]];\n");
    ksprintf(s, "    var lengthMenu10p = [[-1, 10],[\"All\", 10]];\n");
    ksprintf(s, "    var lengthMenu25p = [[-1, 10, 25],[\"All\", 10, 25]];\n");
    ksprintf(s, "    var lengthMenu50p = [[-1, 10, 25, 50],[\"All\", 10, 25, 50]];\n");
    ksprintf(s, "    var lengthMenu100 = [[-1, 10, 25, 50, 100],[\"All\", 10, 25, 50, 100]];\n");
    ksprintf(s, "    var lengthMenu;\n");
    ksprintf(s, "    if(rowCount < 10){\n");
    ksprintf(s, "        lengthMenu = lengthMenuAll;\n");
    ksprintf(s, "    }else if(rowCount < 25){\n");
    ksprintf(s, "        lengthMenu = lengthMenu10p;\n");
    ksprintf(s, "    }else if(rowCount < 50){\n");
    ksprintf(s, "        lengthMenu = lengthMenu25p;\n");
    ksprintf(s, "    }else if(rowCount < 100){\n");
    ksprintf(s, "        lengthMenu = lengthMenu50p;\n");
    ksprintf(s, "    }else{\n");
    ksprintf(s, "        lengthMenu = lengthMenu100;\n");
    ksprintf(s, "    }\n");
    
    ksprintf(s, "    $(document).ready(function() {\n");
    ksprintf(s, "    $('#%s').DataTable( {\n", tid.c_str());
    ksprintf(s, "        lengthMenu: lengthMenu,\n");
    ksprintf(s, "        dom: 'Blfrtip',\n");
    ksprintf(s, "        buttons: [\n");
    ksprintf(s, "            'copyHtml5',\n");
    ksprintf(s, "            {\n");
    ksprintf(s, "                extend: 'excelHtml5',\n");
    ksprintf(s, "                filename: \"%s\",\n", fname.c_str());
    ksprintf(s, "                messageTop: null,\n");
    ksprintf(s, "                messageBottom: null,\n");
    ksprintf(s, "                title: null\n");
    ksprintf(s, "            },\n");
    ksprintf(s, "            {\n");
    ksprintf(s, "                extend: 'csvHtml5',\n");
    ksprintf(s, "                filename: \"%s\",\n", fname.c_str());
    ksprintf(s, "                messageTop: null,\n");
    ksprintf(s, "                messageBottom: null,\n");
    ksprintf(s, "                title: null\n");
    ksprintf(s, "            },\n");
    ksprintf(s, "        ],\n");
    ksprintf(s, "        order: []\n");
    ksprintf(s, "    } );\n");
    ksprintf(s, "    } );\n");
    ksprintf(s, "</script>\n");
}

void HtmlOpt::html2head(kstring_t* s, const std::string& title){
    // beg <html>
    ksprintf(s, "<html>\n");
    // beg <head>
    ksprintf(s, "<head>\n");
    // beg <meta>
    ksprintf(s, "<meta http-equiv=\"content-type\" content=\"text/html;charset=utf-8\" />\n");
    // end <meta>
    // beg <title>
    ksprintf(s, "<title>%s report at %s </title>\n", title.c_str(), util::currentTime().c_str());
    // end <title>
    printJsAndCSSRes(s);
    // beg <script>
    ksprintf(s, "<script type=\"text/javascript\">\n");
    ksprintf(s, "    function showOrHide(divname) {\n");
    ksprintf(s, "        div = document.getElementById(divname);\n");
    ksprintf(s, "        if(div.style.display == 'none')\n");
    ksprintf(s, "            div.style.display = 'block';\n");
    ksprintf(s, "        else\n");
    ksprintf(s, "            div.style.display = 'none';\n");
    ksprintf(s, "    }\n");
    ksprintf(s, "    function toggle(targetid){\n");
    ksprintf(s, "       if (document.getElementById){\n");
    ksprintf(s, "           target=document.getElementById(targetid)\n");
    ksprintf(s, "               if (target.style.display=='table-row'){\n");
    ksprintf(s, "                   target.style.display='none';\n");
    ksprintf(s, "               } else {\n");
    ksprintf(s, "                   target.style.display='table-row';\n");
    ksprintf(s, "               }\n");
    ksprintf(s, "       }\n");
    ksprintf(s, "   }\n");
    ksprintf(s, "   function toggle_target_list(targetid){\n");
    ksprintf(s, "       if (document.getElementById){\n");
    ksprintf(s, "           target=document.getElementById(targetid);\n");
    ksprintf(s, "               if (target.style.display=='block'){\n");
    ksprintf(s, "                   target.style.display='none';\n");
    ksprintf(s, "                   document.getElementById('target_view_btn').value='view';\n");
    ksprintf(s, "               } else {\n");
    ksprintf(s, "                   document.getElementById('target_view_btn').value='hide';\n");
    ksprintf(s, "                   target.style.display='block';\n");
    ksprintf(s, "               }\n");
    ksprintf(s, "       }\n");
    ksprintf(s, "   }\n");
    ksprintf(s, "</script>\n");
    // end <script>
    // beg <style>
    ksprintf(s, "<style type=\"text/css\">\n");
    ksprintf(s, "td {border:1px solid #dddddd;padding:5px;font-size:12px;}\n");
    ksprintf(s, "table {border:1px solid #999999;padding:2x;border-collapse:collapse; width:1000px}\n");
    ksprintf(s, ".col1 {width:300px; font-weight:bold;}\n");
    ksprintf(s, "img {padding:30px;}\n");
    ksprintf(s, "#menu {font-family:Consolas, 'Liberation Mono', Menlo, Courier, monospace;}\n");
    ksprintf(s, "#menu a {color:#0366d6; font-size:18px;font-weight:600;line-height:28px;text-decoration:none;");
    ksprintf(s, "font-family:-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol'}\n");
    ksprintf(s, "a:visited {color: #999999}\n");
    ksprintf(s, ".alignleft {text-align:left;}\n");
    ksprintf(s, ".alignright {text-align:right;}\n");
    ksprintf(s, ".figure {width:%dpx;height:%dpx;}\n", figw, figh);
    ksprintf(s, ".header {color:#ffffff;padding:1px;height:20px;background:#000000;}\n");
    ksprintf(s, ".section_title {color:#ffffff;font-size:20px;padding:5px;text-align:left;background:#000099; margin-top:10px;}\n");
    ksprintf(s, ".section_title_lite {color:#ffffff;font-size:16px;padding:5px;text-align:left;background:#000099; margin-left:10px; margin-top:10px;}\n");
    ksprintf(s, ".section_title_offset_left {color:#ffffff;font-size:16px;padding:5px;text-align:left;background:#000099;margin-left:20px; margin-top:10px;}\n");
    ksprintf(s, ".subsection_title {font-size:16px;padding:5px;margin-top:10px;text-align:left;color:#000099}\n");
    ksprintf(s, "#container {text-align:center;padding:3px 3px 3px 10px;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}\n");
    ksprintf(s, ".menu_item {text-align:left;padding-top:5px;font-size:18px;}\n");
    ksprintf(s, ".highlight {text-align:left;padding-top:30px;padding-bottom:30px;font-size:20px;line-height:35px;}\n");
    ksprintf(s, "#helper {text-align:left;border:1px dotted #fafafa;color:#777777;font-size:12px;}\n");
    ksprintf(s, "#footer {text-align:left;padding:15px;color:#ffffff;font-size:10px;background:#000099;");
    ksprintf(s, "font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}\n");
    ksprintf(s, ".sub_section_tips {color:#999999;font-size:10px;padding-left:5px;padding-bottom:3px;}\n");
    ksprintf(s, "</style>\n");
    // end <style>
    ksprintf(s, "</head>\n");
    // end <head>
    // beg <body>
    ksprintf(s, "<body><div style=\"text-align:left;\" id='container'>\n");
}

void HtmlOpt::html2foot(kstring_t* s){
    // end <body><container...
    ksprintf(s, "</div>\n</body>\n");
    // begin footer
    ksprintf(s, "<div id='footer'> ");
    ksprintf(s, "<p>analysis finished @ %s </p>", util::currentTime().c_str());
    ksprintf(s, "</div>\n");
    // end footer
    // end html
    ksprintf(s, "</html>\n");
}

void HtmlOpt::html2endbody(kstring_t* s){
    // end <body><container...
    ksprintf(s, "</div>\n</body>\n");
}

void HtmlOpt::html2footonly(kstring_t* s){
    // begin footer
    ksprintf(s, "<div id='footer'> ");
    ksprintf(s, "<p>analysis finished @ %s </p>", util::currentTime().c_str());
    ksprintf(s, "</div>\n");
    // end footer
    // end html
    ksprintf(s, "</html>\n");
}

void HtmlOpt::printJssave2svg(kstring_t* s){
    ksprintf(s, "<script type=\"text/javascript\">\n");
    ksprintf(s, "    function expsvg_getSvgElement(svg) {\n");
    ksprintf(s, "        var div = document.createElement(\"div\");\n");
    ksprintf(s, "        div.className = \"tempdiv-svg-exportJS\";\n");
    ksprintf(s, "        if (typeof svg === \"string\") {\n");
    ksprintf(s, "            div.insertAdjacentHTML(\"beforeend\", svg.trim());\n");
    ksprintf(s, "            svg = div.firstChild;\n");
    ksprintf(s, "        } \n");
    ksprintf(s, "        if (!svg.nodeType || svg.nodeType !== 1) {\n");
    ksprintf(s, "            return null;\n");
    ksprintf(s, "        } \n");
    ksprintf(s, "        var svgClone = svg.cloneNode(true);\n");
    ksprintf(s, "        svgClone.style.display = null;\n");
    ksprintf(s, "        div.appendChild(svgClone);\n");
    ksprintf(s, "        div.style.visibility = \"hidden\";\n");
    ksprintf(s, "        div.style.display = \"table\";\n");
    ksprintf(s, "        div.style.position = \"absolute\";\n");
    ksprintf(s, "        document.body.appendChild(div);\n");
    ksprintf(s, "        return svgClone;\n");
    ksprintf(s, "    }\n");
    ksprintf(s, "\n");
    ksprintf(s, "    function expsvg_setOptions(svgElement, options) {\n");
    ksprintf(s, "        _options = {\n");
    ksprintf(s, "            originalWidth: 100,\n");
    ksprintf(s, "            originalHeight: 100,\n");
    ksprintf(s, "            width: 100,\n");
    ksprintf(s, "            height: 100, \n");
    ksprintf(s, "            scale: 1,\n");
    ksprintf(s, "            useCSS: true,\n");
    ksprintf(s, "            transparentBackgroundReplace: \"white\",\n");
    ksprintf(s, "            allowCrossOriginImages: false,\n");
    ksprintf(s, "        };\n");
    ksprintf(s, "        _options.originalHeight = svgElement.style.getPropertyValue(\"height\").indexOf(\"%%\") !== -1 \n");
    ksprintf(s, "            || (svgElement.getAttribute(\"height\") && svgElement.getAttribute(\"height\").indexOf(\"%%\") !== -1 )\n");
    ksprintf(s, "            ? svgElement.getBBox().height * _options.scale\n");
    ksprintf(s, "            : svgElement.getBoundingClientRect().height * _options.scale;\n");
    ksprintf(s, "        _options.originalWidth = svgElement.style.getPropertyValue(\"width\").indexOf(\"%%\") !== -1 \n");
    ksprintf(s, "            || (svgElement.getAttribute(\"width\") && svgElement.getAttribute(\"width\").indexOf(\"%%\") !== -1 )\n");
    ksprintf(s, "            ? svgElement.getBBox().width * _options.scale\n");
    ksprintf(s, "            : svgElement.getBoundingClientRect().width * _options.scale;\n");
    ksprintf(s, "        if (options && options.scale && typeof options.scale === \"number\") {\n");
    ksprintf(s, "            _options.scale = options.scale;\n");
    ksprintf(s, "        }\n");
    ksprintf(s, "        if (!options || !options.height) {\n");
    ksprintf(s, "            _options.height = _options.originalHeight * _options.scale;\n");
    ksprintf(s, "        }\n");
    ksprintf(s, "        else if (typeof options.height === \"number\") {\n");
    ksprintf(s, "            _options.height = options.height * _options.scale;\n");
    ksprintf(s, "        }\n");
    ksprintf(s, "        if (!options || !options.width) {\n");
    ksprintf(s, "            _options.width = _options.originalWidth * _options.scale;\n");
    ksprintf(s, "        }\n");
    ksprintf(s, "        else if (typeof options.width === \"number\") {\n");
    ksprintf(s, "            _options.width = options.width * _options.scale;\n");
    ksprintf(s, "        } \n");
    ksprintf(s, "        if (options && options.useCSS === false) {\n");
    ksprintf(s, "            _options.useCSS = false;\n");
    ksprintf(s, "        }\n");
    ksprintf(s, "        if (options && options.transparentBackgroundReplace) {\n");
    ksprintf(s, "            _options.transparentBackgroundReplace = options.transparentBackgroundReplace;\n");
    ksprintf(s, "        }\n");
    ksprintf(s, "        if (options && options.allowCrossOriginImages) {\n");
    ksprintf(s, "            _options.allowCrossOriginImages = options.allowCrossOriginImages;\n");
    ksprintf(s, "        }\n");
    ksprintf(s, "    }\n");
    ksprintf(s, "\n");
    ksprintf(s, "    function expsvg_useCSSfromComputedStyles(element, elementClone) {\n");
    ksprintf(s, "        if (typeof getComputedStyle !== \"function\"){\n");
    ksprintf(s, "            return;\n");
    ksprintf(s, "        } \n");
    ksprintf(s, "        element.childNodes.forEach(function(child, index){\n");
    ksprintf(s, "            if (child.nodeType === 1/*Node.ELEMENT_NODE*/) {\n");
    ksprintf(s, "                expsvg_useCSSfromComputedStyles(child, elementClone.childNodes[parseInt(index, 10)]);\n");
    ksprintf(s, "            }\n");
    ksprintf(s, "        });\n");
    ksprintf(s, "        \n");
    ksprintf(s, "        var compStyles = window.getComputedStyle(element);\n");
    ksprintf(s, "        if (compStyles.length > 0) {\n");
    ksprintf(s, "            for (const compStyle of compStyles){\n");
    ksprintf(s, "                if ([\"width\", \"height\", \"inline-size\", \"block-size\"].indexOf(compStyle) === -1 ) {\n");
    ksprintf(s, "                    elementClone.style.setProperty(compStyle, compStyles.getPropertyValue(compStyle));\n");
    ksprintf(s, "                }\n");
    ksprintf(s, "            };\n");
    ksprintf(s, "        }\n");
    ksprintf(s, "    }\n");
    ksprintf(s, "\n");
    ksprintf(s, "    function expsvg_setupSvg(svgElement, originalSvg, asString)\n");
    ksprintf(s, "    {\n");
    ksprintf(s, "        if (typeof asString === \"undefined\") { asString = true; }\n");
    ksprintf(s, "        if (_options.useCSS && typeof originalSvg === \"object\") {\n");
    ksprintf(s, "            expsvg_useCSSfromComputedStyles(originalSvg, svgElement);\n");
    ksprintf(s, "            svgElement.style.display = null;\n");
    ksprintf(s, "        }\n");
    ksprintf(s, "\n");
    ksprintf(s, "        svgElement.style.width = _options.width;\n");
    ksprintf(s, "        svgElement.style.height = _options.height;\n");
    ksprintf(s, "        svgElement.setAttribute(\"width\", _options.width);\n");
    ksprintf(s, "        svgElement.setAttribute(\"height\", _options.height);\n");
    ksprintf(s, "        svgElement.setAttribute(\"preserveAspectRatio\", \"none\");\n");
    ksprintf(s, "        svgElement.setAttribute(\"viewBox\", \"0 0 \" + (_options.originalWidth) + \" \" + (_options.originalHeight));\n");
    ksprintf(s, "       \n");
    ksprintf(s, "        var elements = document.getElementsByClassName(\"tempdiv-svg-exportJS\");\n");
    ksprintf(s, "        while(elements.length > 0){\n");
    ksprintf(s, "            elements[0].parentNode.removeChild(elements[0]);\n");
    ksprintf(s, "        }\n");
    ksprintf(s, "        //get svg string\n");
    ksprintf(s, "        if (asString)\n");
    ksprintf(s, "        {\n");
    ksprintf(s, "            var serializer = new XMLSerializer();\n");
    ksprintf(s, "            //setting currentColor to black matters if computed styles are not used\n");
    ksprintf(s, "            var svgString = serializer.serializeToString(svgElement).replace(/currentColor/g, \"black\");\n");
    ksprintf(s, "\n");
    ksprintf(s, "            //add namespaces\n");
    ksprintf(s, "            if (!svgString.match(/^<svg[^>]+xmlns=\"http\\:\\/\\/www\\.w3\\.org\\/2000\\/svg\"/)) {\n");
    ksprintf(s, "                svgString = svgString.replace(/^<svg/, \"<svg xmlns=\\\"http://www.w3.org/2000/svg\\\"\");\n");
    ksprintf(s, "            }\n");
    ksprintf(s, "            if (!svgString.match(/^<svg[^>]+\"http\\:\\/\\/www\\.w3\\.org\\/1999\\/xlink\"/)) {\n");
    ksprintf(s, "                svgString = svgString.replace(/^<svg/, \"<svg xmlns:xlink=\\\"http://www.w3.org/1999/xlink\\\"\");\n");
    ksprintf(s, "            }\n");
    ksprintf(s, "    \n");
    ksprintf(s, "            return svgString;\n");
    ksprintf(s, "        }\n");
    ksprintf(s, "        return svgElement;\n");
    ksprintf(s, "    }\n");
    ksprintf(s, "\n");
    ksprintf(s, "    function expsvg_convertImageURLtoDataURI(image) {\n");
    ksprintf(s, "        return new Promise(function(resolve, reject) {\n");
    ksprintf(s, "            var newImage = new Image();\n");
    ksprintf(s, "                        \n");
    ksprintf(s, "            newImage.onload = function () {\n");
    ksprintf(s, "                var canvas = document.createElement(\"canvas\");\n");
    ksprintf(s, "                canvas.width = this.naturalWidth || this.getAttribute(\"width\") || this.style.getPropertyValue(\"width\") || 300; \n");
    ksprintf(s, "                canvas.height = this.naturalHeight || this.getAttribute(\"height\") || this.style.getPropertyValue(\"height\") || 300; \n");
    ksprintf(s, "\n");
    ksprintf(s, "                canvas.getContext(\"2d\").drawImage(this, 0, 0);\n");
    ksprintf(s, "\n");
    ksprintf(s, "                var dataURI = canvas.toDataURL(\"image/png\");\n");
    ksprintf(s, "                image.setAttribute(\"href\", dataURI);\n");
    ksprintf(s, "                resolve();\n");
    ksprintf(s, "            };\n");
    ksprintf(s, "            if (_options.allowCrossOriginImages)\n");
    ksprintf(s, "                newImage.crossOrigin = \"anonymous\";\n");
    ksprintf(s, "            newImage.src = image.getAttribute(\"href\") || image.getAttributeNS(\"http://www.w3.org/1999/xlink\", \"href\");\n");
    ksprintf(s, "        });\n");
    ksprintf(s, "    }\n");
    ksprintf(s, "\n");
    ksprintf(s, "    function expsvg_triggerDownload(uri, name, canvas) {\n");
    ksprintf(s, "        name = name.replace(/[/\\?%%*:|\"<>]/g, \"_\");\n");
    ksprintf(s, "        if (navigator.msSaveBlob) {\n");
    ksprintf(s, "            var binary = (decodeURIComponent(uri.split(\",\")[1])), array = [];\n");
    ksprintf(s, "            var mimeString = uri.split(\",\")[0].split(\":\")[1].split(\";\")[0];\n");
    ksprintf(s, "            for (var i = 0; i < binary.length; i++) {\n");
    ksprintf(s, "                array.push(binary.charCodeAt(i));\n");
    ksprintf(s, "            }\n");
    ksprintf(s, "            var blob = null;\n");
    ksprintf(s, "            if (canvas != null) {\n");
    ksprintf(s, "                blob = canvas.msToBlob();\n");
    ksprintf(s, "            } else {\n");
    ksprintf(s, "                blob = new Blob([new Uint8Array(array)], { type: mimeString });\n");
    ksprintf(s, "            }\n");
    ksprintf(s, "            return navigator.msSaveBlob(blob, name);\n");
    ksprintf(s, "        } else {\n");
    ksprintf(s, "            var link = document.createElement(\"a\");\n");
    ksprintf(s, "            link.download = name;\n");
    ksprintf(s, "            link.href = uri;\n");
    ksprintf(s, "            document.body.appendChild(link);\n");
    ksprintf(s, "            link.click();\n");
    ksprintf(s, "            document.body.removeChild(link);\n");
    ksprintf(s, "        }\n");
    ksprintf(s, "    }\n");
    ksprintf(s, "\n");
    ksprintf(s, "    function expsvg_downloadSvg(svg, svgName, options) {\n");
    ksprintf(s, "        var svgElement = expsvg_getSvgElement(svg);\n");
    ksprintf(s, "        if (!svgElement) { return; }\n");
    ksprintf(s, "        if (svgName == null) {\n");
    ksprintf(s, "            svgName = \"chart\";\n");
    ksprintf(s, "        }\n");
    ksprintf(s, "        expsvg_setOptions(svgElement, options);\n");
    ksprintf(s, "        var images = svgElement.getElementsByTagName(\"image\");\n");
    ksprintf(s, "        var image_promises = [];\n");
    ksprintf(s, "        if (images){\n");
    ksprintf(s, "            for (var image of images) {\n");
    ksprintf(s, "                if ((image.getAttribute(\"href\") && image.getAttribute(\"href\").indexOf(\"data:\") === -1)\n");
    ksprintf(s, "                || (image.getAttribute(\"xlink:href\") && image.getAttribute(\"xlink:href\").indexOf(\"data:\") === -1)) {\n");
    ksprintf(s, "                    image_promises.push(expsvg_convertImageURLtoDataURI(image));\n");
    ksprintf(s, "                }\n");
    ksprintf(s, "            }\n");
    ksprintf(s, "        }\n");
    ksprintf(s, "        Promise.all(image_promises).then(function() {\n");
    ksprintf(s, "            var svgString = expsvg_setupSvg(svgElement, svg);\n");
    ksprintf(s, "            svgString = \"<?xml version=\\\"1.0\\\" standalone=\\\"no\\\"?>\\\r\\\n\" + svgString;\n");
    ksprintf(s, "            var url = \"data:image/svg+xml;charset=utf-8,\" + encodeURIComponent(svgString);\n");
    ksprintf(s, "            expsvg_triggerDownload(url, svgName + \".svg\");\n");
    ksprintf(s, "        });\n");
    ksprintf(s, "    }\n");
    ksprintf(s, "    function save2svg(svg_img_id, svg_file_name, svg_width, svg_height){\n");
    ksprintf(s, "        expsvg_downloadSvg(\n");
    ksprintf(s, "        document.getElementById(svg_img_id),\n");
    ksprintf(s, "        svg_file_name,\n");
    ksprintf(s, "        { width: svg_width, height: svg_height}\n");
    ksprintf(s, "        );\n");
    ksprintf(s, "    }\n");
    ksprintf(s, "</script>\n");
}
