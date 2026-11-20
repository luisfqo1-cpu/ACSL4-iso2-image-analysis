import numpy as np

from imagetools_napari import ExampleQWidget


# make_napari_viewer is a pytest fixture that returns a napari viewer object
def test_threshold_slider_adds_thresholded_layer(make_napari_viewer):
    viewer = make_napari_viewer()
    viewer.add_image(np.random.randint(0, 1000, (100, 100)))

    widget = ExampleQWidget(viewer)
    widget.slider.setValue(500)

    assert "thresholded" in [layer.name for layer in viewer.layers]
    np.testing.assert_array_equal(
        viewer.layers["thresholded"].data,
        viewer.layers[0].data > widget.slider.value(),
    )
