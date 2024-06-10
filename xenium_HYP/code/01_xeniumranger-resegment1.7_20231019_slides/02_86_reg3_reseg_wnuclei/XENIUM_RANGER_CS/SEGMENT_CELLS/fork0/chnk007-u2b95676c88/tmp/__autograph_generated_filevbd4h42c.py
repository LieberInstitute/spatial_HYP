# coding=utf-8
from __future__ import annotations
def outer_factory():

    def inner_factory(ag__):

        def tf___infer(self, image: tf.Tensor, offset: tf.Tensor, norm: tf.Tensor, target_height: tf.Tensor, target_width: tf.Tensor) -> tuple[tf.Tensor, tf.Tensor]:
            """Model inference.

        Args:
            image (tf.Tensor): Image of shape (1, height, width, channel). Do not support batch of images at
                the moment.
            offset (tf.Tensor): offset used in normalization
            norm (tf.Tensor): denominator used in normalization
            target_height (tf.Tensor): targeted height for inference. Scale to it if it's different
                with the original height
            target_width (tf.Tensor): targeted width for inference. Scale to it if it's different
                with the original width

        Returns:
            tf.Tensor: row probability of shape (target_height, target_width, 4)
            tf.Tensor: col probability of shape (target_height, target_width, 4)
        """
            with ag__.FunctionScope('_infer', 'fscope', ag__.ConversionOptions(recursive=True, user_requested=True, optional_features=(), internal_convert_user_code=True)) as fscope:
                do_return = False
                retval_ = ag__.UndefinedReturnValue()
                ag__.converted_call(ag__.ld(printwhen), ('Segmentation inference tracing',), None, fscope)
                image = ag__.converted_call(ag__.ld(tf).cond, (ag__.converted_call(ag__.ld(tf).logical_and, (ag__.converted_call(ag__.ld(tf).equal, (ag__.converted_call(ag__.ld(tf).shape, (ag__.ld(image),), None, fscope)[1], ag__.ld(target_height)), None, fscope), ag__.converted_call(ag__.ld(tf).equal, (ag__.converted_call(ag__.ld(tf).shape, (ag__.ld(image),), None, fscope)[2], ag__.ld(target_width)), None, fscope)), None, fscope),), dict(true_fn=ag__.autograph_artifact(lambda : ag__.ld(image)), false_fn=ag__.autograph_artifact(lambda : ag__.converted_call(ag__.ld(tf).image.resize, (ag__.ld(image), [ag__.ld(target_height), ag__.ld(target_width)]), dict(method=ag__.ld(tf).image.ResizeMethod.BILINEAR, preserve_aspect_ratio=False), fscope))), fscope)
                image = (ag__.ld(image) - ag__.ld(offset)) / ag__.ld(norm)
                output = ag__.converted_call(ag__.ld(tf).graph_util.import_graph_def, (ag__.ld(self)._graph_def,), dict(name='', input_map={'args_0': ag__.ld(image)}, return_elements=['Identity:0']), fscope)[0]
                (direct_map, stationary_px) = ag__.converted_call(ag__.ld(tf).graph_util.import_graph_def, (ag__.ld(self)._post_model,), dict(name='', input_map={'model_output': ag__.ld(output)}, return_elements=['Identity:0', 'Identity_1:0']), fscope)
                try:
                    do_return = True
                    retval_ = (ag__.ld(direct_map), ag__.ld(stationary_px))
                except:
                    do_return = False
                    raise
                return fscope.ret(retval_, do_return)
        return tf___infer
    return inner_factory